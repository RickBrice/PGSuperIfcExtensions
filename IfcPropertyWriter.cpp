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

   void WritePropertySets(ElementKind role, const std::vector<IfcExport::Schema::IfcObjectDefinition>& objects, const ExportContext& context) override
   {
      ::WritePropertySets<IfcExport::Schema>(Hooks().file, role, objects, WithSession(context));
   }

   std::vector<IfcExport::Schema::IfcPropertySetDefinition> CreateTypePropertySets(ElementKind role, const ExportContext& context) override
   {
      return ::CreateTypePropertySets<IfcExport::Schema>(Hooks().file, role, WithSession(context));
   }

   void Classify(ElementKind role, IfcExport::Schema::IfcObjectDefinition object) override
   {
      ::Classify<IfcExport::Schema>(Hooks().file, role, object, Hooks().options);
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

   // the context with the broker and options of the export, if the extension left them out
   ExportContext WithSession(const ExportContext& context) const
   {
      ExportContext c(context);
      if (!c.broker) c.broker = Hooks().broker;
      if (!c.options) c.options = &Hooks().options;
      return c;
   }
};

CIfcExportSession::CIfcExportSession(const CIfcExportOptions& options, IIfcExportExtension* extension)
{
   ASSERT(ms_pCurrent == nullptr); // one export at a time

   if (extension)
      m_pExtensionContent = RegisterExtension(*extension);
   m_pScope = std::make_unique<CIfcExtensionScope>(m_pExtensionContent.get());

   std::filesystem::path path(options.mapping_file.GetString());
   m_pTable = CIfcMappingTable::LoadActive(path);

   if (extension)
      m_pHooks = extension->CreateHooks();

   ms_pCurrent = this;
}

CIfcExportSession::~CIfcExportSession()
{
   if (m_pHooks)
      m_pHooks->OnExportEnd(m_bSucceeded);

   ms_pCurrent = nullptr;
}

void CIfcExportSession::StartHooks(IfcExport::File& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options)
{
   if (!m_pHooks)
      return;

   m_pServices = std::make_unique<CIfcExportServices>(*this);
   m_pHookContext = std::make_unique<IfcExport::HookContext>(IfcExport::HookContext{ file, pBroker, options, *m_pServices, m_Elements });
}

void CIfcExportSession::CallHook(const char* hook, const std::function<void(IfcExport::IHooks&, IfcExport::HookContext&)>& call)
{
   if (!m_pHookContext)
      return;

   auto failed = [this, hook](const std::string& problem)
      {
         return CIfcExtensionException("The export extension \"" + m_pExtensionContent->id + "\" failed in " + hook + ": " + problem);
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
   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   switch (target.display_unit)
   {
   case ExportUnit::SpanLength: value = WBFL::Units::ConvertFromSysUnits(value, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure); break;
   case ExportUnit::Deflection: value = WBFL::Units::ConvertFromSysUnits(value, pDisplayUnits->GetDeflectionUnit().UnitOfMeasure); break;
   case ExportUnit::Stress: value = WBFL::Units::ConvertFromSysUnits(value, pDisplayUnits->GetStressUnit().UnitOfMeasure); break;
   case ExportUnit::Angle: value = WBFL::Units::ConvertFromSysUnits(value, pDisplayUnits->GetAngleUnit().UnitOfMeasure); break;
   case ExportUnit::SmallArea: value = WBFL::Units::ConvertFromSysUnits(value, pDisplayUnits->GetAreaUnit().UnitOfMeasure); break;
   case ExportUnit::BigArea: value = WBFL::Units::ConvertFromSysUnits(value, WBFL::Units::Measure::Feet2); break; // as GetBigAreaUnit
   case ExportUnit::Mass: value = WBFL::Units::ConvertFromSysUnits(value, WBFL::Units::Measure::PoundMass); break; // as GetMassUnit
   default: break;
   }

   if (0 < target.display_round)
      value = RoundOff(value, target.display_round);

   return value;
}
