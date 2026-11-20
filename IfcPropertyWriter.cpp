///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright © 1999-2026  Washington State Department of Transportation
//                        Bridge and Structures Office
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the Alternate Route Open Source License as
// published by the Washington State Department of Transportation,
// Bridge and Structures Office.
//
// This program is distributed in the hope that it will be useful, but
// distribution is AS IS, WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
// the Alternate Route Open Source License for more details.
//
// You should have received a copy of the Alternate Route Open Source
// License along with this program; if not, write to the Washington
// State Department of Transportation, Bridge and Structures Office,
// P.O. Box  47340, Olympia, WA 98503, USA or e-mail
// Bridge_Support@wsdot.wa.gov
///////////////////////////////////////////////////////////////////////
#include "stdafx.h"
#include "IfcPropertyWriter.h"

#include <sstream>

CIfcExportSession* CIfcExportSession::ms_pCurrent = nullptr;

// What an export extension's hooks use to write property sets and log (IfcExport::IServices). The file, broker, and options
// are the hook context's
class CIfcExportServices : public IfcExport::IServices
{
public:
   explicit CIfcExportServices(CIfcExportSession& session) : m_Session(session)
   {
   }

   ElementKind Role(const std::string& name) const override
   {
      ElementKind kind;
      if (!GetElementKind(name, kind))
         throw CIfcExtensionException("unknown element role \"" + name + "\"");
      return kind;
   }

   bool HasPropertySets(ElementKind role) const override
   {
      const auto& table = m_Session.GetTable();
      return !table.GetPropertySets(role, PropertyOwner::Occurrence).empty() || !table.GetPropertySets(role, PropertyOwner::Type).empty() || !table.GetClassifications(role).empty();
   }

   void WritePropertySets(ElementKind role, const std::vector<IfcExport::Schema::IfcObjectDefinition>& objects, const IfcExport::ElementKeys& keys) override
   {
      ::WritePropertySets<IfcExport::Schema>(Hooks().file, role, objects, Context(keys));
   }

   std::vector<IfcExport::Schema::IfcPropertySetDefinition> CreateTypePropertySets(ElementKind role, const IfcExport::ElementKeys& keys) override
   {
      return ::CreateTypePropertySets<IfcExport::Schema>(Hooks().file, role, Context(keys));
   }

   IfcExport::Schema::IfcMaterial RebarMaterial(const WBFL::Materials::Rebar& rebar) override
   {
      return ::GetExportRebarMaterial(Hooks().file, Hooks().broker, Options(), rebar);
   }

   void Classify(ElementKind role, IfcExport::Schema::IfcObjectDefinition object) override
   {
      ::Classify<IfcExport::Schema>(Hooks().file, role, object, Options());
   }

   IfcExport::Schema::IfcGeometricRepresentationSubContext BodySubContext() override
   {
      return Hooks().file.getRepresentationSubContext(std::string("Body"), std::string("Model"));
   }

   void Log(const std::string& line) override
   {
      m_Session.m_HookLog.push_back(line);
   }

private:
   CIfcExportSession& m_Session;

   const IfcExport::HookContext& Hooks() const { return *m_Session.m_pHookContext; }
   const CIfcExportOptions& Options() const { return *m_Session.m_pOptions; }

   // the export context of an extension's element: the export's broker and options, and the element's keys
   ExportContext Context(const IfcExport::ElementKeys& keys) const
   {
      ExportContext c;
      c.broker = Hooks().broker;
      c.options = &Options();
      c.segment = keys.segment;
      c.pier = keys.pier;
      c.span = keys.span;
      c.index = keys.index;
      c.extension = keys.extension;
      return c;
   }
};

CIfcExportSession::CIfcExportSession(const CIfcExportOptions& options, const CExportProfile& profile) :
   m_Profile(profile.Describe())
{
   ASSERT(ms_pCurrent == nullptr); // one export at a time

   if (auto* extension = profile.GetExtension())
      m_pExtensionContent = RegisterExtension(*extension);
   m_pScope = std::make_unique<CIfcExtensionScope>(m_pExtensionContent.get());

   std::filesystem::path path(options.mapping_file.GetString());
   m_pTable = CIfcMappingTable::LoadActive(path, options.mapping_source);

   m_pHooks = profile.CreateHooks();

   ms_pCurrent = this;
}

CIfcExportSession::~CIfcExportSession()
{
   m_pHooks->OnExportEnd(m_bSucceeded);
   ms_pCurrent = nullptr;
}

void CIfcExportSession::StartHooks(IfcExport::File& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options)
{
   m_pOptions = &options;
   m_BuildOptions.mapping_file = options.mapping_file;
   m_BuildOptions.include_rebar = options.include_rebar;
   m_BuildOptions.include_quantities = options.include_quantities;
   m_BuildOptions.display_units_for_properties = options.display_units_for_properties;

   m_pServices = std::make_unique<CIfcExportServices>(*this);
   m_pHookContext = std::make_unique<IfcExport::HookContext>(IfcExport::HookContext{ file, pBroker, m_BuildOptions, *m_pServices, m_Elements });
}

void CIfcExportSession::CallHook(const char* hook, const std::function<void(IfcExport::IHooks&, IfcExport::HookContext&)>& call)
{
   if (!m_pHookContext)
      return;

   auto failed = [this, hook](const std::string& problem)
      {
         // e.g. "The export extension "MyAgency.MyExchanger" failed in OnPier: ..."
         std::string who = m_Profile;
         if (!who.empty())
            who[0] = (char)toupper((unsigned char)who[0]);
         return CIfcExtensionException(who + " failed in " + hook + ": " + problem);
      };

   try
   {
      call(*m_pHooks, *m_pHookContext);
   }
   catch (const std::exception& e) // CIfcExtensionException too
   {
      throw failed(e.what());
   }
   catch (...)
   {
      throw failed("unexpected error");
   }
}

CIfcExportSession& CIfcExportSession::Current()
{
   ASSERT(ms_pCurrent); // only available during an export
   return *ms_pCurrent;
}

bool IncludePropertySet(const PropertySetDeclaration& pset, const CIfcExportOptions& options)
{
   switch (pset.condition)
   {
   case PropertySetDeclaration::Condition::Classify: return options.classify;
   case PropertySetDeclaration::Condition::Quantities: return options.include_quantities;
   case PropertySetDeclaration::Condition::Always: return true;
   }
   ASSERT(false);
   return true;
}

bool UseDisplayUnits(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options)
{
   if (!options.display_units_for_properties)
      return false;

   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   return pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US;
}

Float64 ConvertToDisplayUnits(std::shared_ptr<WBFL::EAF::Broker> pBroker, const TargetDef& target, Float64 value)
{
   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);
   switch (target.display_unit)
   {
   case ExportUnit::SpanLength: value = WBFL::Units::ConvertFromSysUnits(value, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure); break;
   case ExportUnit::Deflection: value = WBFL::Units::ConvertFromSysUnits(value, pDisplayUnits->GetDeflectionUnit().UnitOfMeasure); break;
   case ExportUnit::Stress: value = WBFL::Units::ConvertFromSysUnits(value, pDisplayUnits->GetStressUnit().UnitOfMeasure); break;
   case ExportUnit::Angle: value = WBFL::Units::ConvertFromSysUnits(value, pDisplayUnits->GetAngleUnit().UnitOfMeasure); break;
   case ExportUnit::SmallArea: value = WBFL::Units::ConvertFromSysUnits(value, pDisplayUnits->GetAreaUnit().UnitOfMeasure); break;
   case ExportUnit::BigArea: value = WBFL::Units::ConvertFromSysUnits(value, WBFL::Units::Measure::Feet2); break; // as GetBigAreaUnit
   case ExportUnit::Mass: value = WBFL::Units::ConvertFromSysUnits(value, WBFL::Units::Measure::PoundMass); break; // as GetMassUnit
   case ExportUnit::Volume: value = WBFL::Units::ConvertFromSysUnits(value, WBFL::Units::Measure::Feet3); break; // as GetVolumeUnit
   case ExportUnit::ComponentDim: value = WBFL::Units::ConvertFromSysUnits(value, pDisplayUnits->GetComponentDimUnit().UnitOfMeasure); break;
   default: break;
   }

   if (0 < target.display_round)
      value = RoundOff(value, target.display_round);

   return value;
}

ExportValue TargetExportValue(const TargetDef& target, const ExportContext& context)
{
   if (target.needs_analysis && !context.options->include_analysis_values)
      return ExportValue::NoValue();

   return target.get(context);
}

std::string GetTypePropertySignature(ElementKind role, const ExportContext& context)
{
   std::ostringstream os;
   for (const auto* pset : CIfcExportSession::Current().GetTable().GetPropertySets(role, PropertyOwner::Type))
   {
      if (!IncludePropertySet(*pset, *context.options))
         continue;

      for (const auto& property : pset->properties)
      {
         if (property.value || !property.target || !property.target->get)
            continue; // the same for every element

         auto result = TargetExportValue(*property.target, context);
         os << pset->name << "." << property.name << "=" << (int)result.state;
         if (result.state == ExportValue::State::Value)
         {
            if (result.list.empty())
               os << ":" << FormatTargetValue(result.value);
            for (const auto& v : result.list)
               os << ":" << FormatTargetValue(v);
         }
         os << ";";
      }
   }
   return os.str();
}
