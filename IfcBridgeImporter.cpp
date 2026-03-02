///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright © 1999-2025  Washington State Department of Transportation
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
#include "IfcBridgeImporter.h"
#include "IfcImporter.h"
#include "Properties.h"
#include "IfcImporterException.h"
#include "USBridge_Classifications.h"
#include "Utilities.h"

#include <MFCTools\Prompts.h>
#include <boost/range/combine.hpp>
#include <psgLib/BridgeDescription2.h>
#include <psgLib/GirderLabel.h>
#include <psgLib/GirderLibraryEntry.h>


std::pair<GroupIndexType, GirderIndexType> ExtractSpanAndGirder(const std::string& s)
{
   GroupIndexType grpIdx = INVALID_INDEX;
   GirderIndexType gdrIdx = INVALID_INDEX;
   std::istringstream iss(s);
   std::string word;

   while (iss >> word)
   {
      if (word == "Span")
         iss >> grpIdx;
      else if (word == "Girder")
         iss >> gdrIdx;
      else if (word == ",")
      { // do nothing
      }
      else
      {
         USES_CONVERSION;
         std::_tostringstream os;
         os << _T("Unexpected DesignLocationNumber property in Pset_PrecastConcreteElementGeneral property set (") << A2T(s.c_str()) << _T(")");
         IFC_THROW(os.str().c_str());
      }
   }

   return { grpIdx - 1,gdrIdx - 1 };
}

Ifc4x3_add2::IfcBridgePart* GetBridgePart(IfcParse::IfcFile& file, Ifc4x3_add2::IfcBridgePartTypeEnum part_type)
{
   auto parts = file.instances_by_type<Ifc4x3_add2::IfcBridgePart>();
   for (auto part : *parts)
   {
      if (part->PredefinedType().has_value() && part->PredefinedType().get() == part_type)
      {
         return part;
      }
   }
   return nullptr;
}

std::vector<Ifc4x3_add2::IfcBridgePart*> GetBridgeParts(IfcParse::IfcFile& file, Ifc4x3_add2::IfcBridgePartTypeEnum part_type)
{
   std::vector<Ifc4x3_add2::IfcBridgePart*> parts_found;
   auto parts = file.instances_by_type<Ifc4x3_add2::IfcBridgePart>();
   for (auto part : *parts)
   {
      if (part->PredefinedType().has_value() && part->PredefinedType().get() == part_type)
      {
         parts_found.push_back(part);
      }
   }
   return parts_found;
}

template <typename Schema>
typename Schema::IfcMaterial* GetMaterial(typename Schema::IfcObjectDefinition* objectdef)
{
   auto associations = objectdef->HasAssociations();
   if (associations)
   {
      for (auto rel : *associations)
      {
         auto rel_associates_material = rel->as<Ifc4x3_add2::IfcRelAssociatesMaterial>();
         if (rel_associates_material)
         {
            auto material = rel_associates_material->RelatingMaterial();
            return material->as<typename Schema::IfcMaterial>();
         }
      }
   }

   return nullptr;
}

int GetBeamTypeCount(IfcParse::IfcFile& file)
{
   int count = 0;
   auto beam_types = file.instances_by_type<Ifc4x3_add2::IfcBeamType>();
   for (auto beam_type : *beam_types)
   {
      if (beam_type->PredefinedType() == Ifc4x3_add2::IfcBeamTypeEnum::IfcBeamType_BEAM)
         count++;
   }

   return count;
}


CIfcBridgeImporter::CIfcBridgeImporter(CIfcImporter& importer) :
   m_Importer(importer)
{
}

CIfcImporter::ImportResult CIfcBridgeImporter::Import(IfcParse::IfcFile& file)
{
   auto bridge = GetBridge(file);
   if (bridge == nullptr)
      return CIfcImporter::ImportResult::NotFound;

   auto parts = file.instances_by_type<Ifc4x3_add2::IfcBridgePart>();
   PierIndexType nPiers = 0;
   for (auto& part : *parts)
   {
      if (part->PredefinedType().has_value() && (part->PredefinedType().get() == Ifc4x3_add2::IfcBridgePartTypeEnum::IfcBridgePartType_ABUTMENT || part->PredefinedType().get() == Ifc4x3_add2::IfcBridgePartTypeEnum::IfcBridgePartType_PIER))
         nPiers++;
   }


   SpanIndexType nSpans = INVALID_INDEX;
   auto value = GetProperty<Ifc4x3_add2, Ifc4x3_add2::IfcInteger>(bridge, "usBridge_BridgeCommon", "usBridge_NumberOfSpans");
   if (value)
   {
      nSpans = (SpanIndexType)(*value);
      if (nSpans != nPiers - 1)
         IFC_THROW(_T("Number of spans modeled does not match number of spans in usBridge_BridgeCommon property set"));
   }
   else
   {
      IFC_THROW(_T("usBridge_NumberOfSpans property not found in usBridge_BridgeCommon property set"));
   }

   std::vector<GirderIndexType> nGirders;
   nGirders.assign(nSpans, 0);
   auto superstructure = GetBridgePart(file, Ifc4x3_add2::IfcBridgePartTypeEnum::Value::IfcBridgePartType_SUPERSTRUCTURE);
   auto rel_contained_elements = superstructure->ContainsElements();
   for (auto contained_element : *rel_contained_elements)
   {
      auto related_elements = contained_element->RelatedElements();
      for (auto related_element : *related_elements)
      {
         auto beam = related_element->as<Ifc4x3_add2::IfcBeam>();
         if (beam)
         {
            auto value = GetProperty<Ifc4x3_add2, Ifc4x3_add2::IfcLabel>(beam, "Pset_PrecastConcreteElementGeneral", "DesignLocationNumber");
            if (!value)
            {
               IFC_THROW(_T("DesignLocationNumber property not found in Pset_PrecastConcreteElementGeneral"));
            }
            auto [spanIdx, gdrIdx] = ExtractSpanAndGirder(*value);
            nGirders[spanIdx]++;
         }
      }
   }

   bool bSameNumGirdersInAllSpans = std::adjacent_find(nGirders.begin(), nGirders.end(), std::not_equal_to<>()) == nGirders.end() ? true : false;

   // Get the existing bridge description. We are going to modify the bridge description with
   // information extracted from the IFC file.
   GET_IFACE2(m_Importer.GetBroker(),IBridgeDescription, pIBridgeDesc);
   auto bridge_desc = *(pIBridgeDesc->GetBridgeDescription());

   SpanIndexType nSpansToAdd = nSpans - pIBridgeDesc->GetSpanCount();
   if (0 < nSpansToAdd)
   {
      for (SpanIndexType i = 0; i < nSpansToAdd; i++)
      {
         bridge_desc.AppendSpan(nullptr, nullptr, true, 0);
      }

      bridge_desc.UseSameGirderForEntireBridge(true);

      if (bSameNumGirdersInAllSpans)
      {
         bridge_desc.UseSameNumberOfGirdersInAllGroups(true);
         bridge_desc.SetGirderCount(nGirders.front());
      }
   }

   if (!bSameNumGirdersInAllSpans)
   {
      bridge_desc.UseSameNumberOfGirdersInAllGroups(false);
      for (SpanIndexType spanIdx = 0; spanIdx < nSpans; spanIdx++)
      {
         bridge_desc.GetGirderGroup(spanIdx)->SetGirderCount(nGirders[spanIdx]);
      }
   }

   //
   // Position the abutments and piers
   //

   // Per TPF modeling guide, piers and abutments are different types.
   // Get the abutments and piers and put into a single vector because we need to treat them the same in PGSuper
   std::vector<Ifc4x3_add2::IfcBridgePart*> abutments = GetBridgeParts(file, Ifc4x3_add2::IfcBridgePartTypeEnum::Value::IfcBridgePartType_ABUTMENT);
   std::vector<Ifc4x3_add2::IfcBridgePart*> piers = GetBridgeParts(file, Ifc4x3_add2::IfcBridgePartTypeEnum::Value::IfcBridgePartType_PIER);
   piers.insert(piers.begin(), abutments.front());
   piers.insert(piers.end(), abutments.back());

   nPiers = bridge_desc.GetPierCount();
   for (PierIndexType pierIdx = 0; pierIdx < nPiers; pierIdx++)
   {
      // get the pier station from the positioning element and set it on the PGSuper pier
      auto pPier = bridge_desc.GetPier(pierIdx);
      auto pier = piers[pierIdx];
      auto rel_positions = pier->PositionedRelativeTo();
      if (rel_positions->size() == 0)
      {
         std::_tostringstream os;
         os << _T("Pier ") << LABEL_PIER(pierIdx) << _T(": must be positioned with an IfcReferent");
         IFC_THROW(os.str().c_str());
      }

      auto positioning_element = (*rel_positions->begin())->RelatingPositioningElement();
      auto ref = positioning_element->as<Ifc4x3_add2::IfcReferent>();
      if (!ref)
      {
         std::_tostringstream os;
         os << _T("Pier ") << LABEL_PIER(pierIdx) << _T(": must be positioned with an IfcReferent");
         IFC_THROW(os.str().c_str());
      }

      auto station = GetProperty<Ifc4x3_add2, Ifc4x3_add2::IfcLengthMeasure>(ref, "Pset_Stationing", "Station");
      if (!station)
      {
         std::_tostringstream os;
         os << _T("Pier ") << LABEL_PIER(pierIdx) << _T(": Station property in Pset_Stationing not found");
         IFC_THROW(os.str().c_str());
      }

      pPier->SetStation(*station);
   }

   // NOTE: This is not the cleanest way to do this, but it gets the job done for now.
   // I want to keep girder spacing from a custom property set separate from the pier stationing.
   bridge_desc.SetGirderSpacingType(pgsTypes::SupportedBeamSpacing::sbsGeneral);
   for (PierIndexType pierIdx = 0; pierIdx < nPiers; pierIdx++)
   {
      auto pier = piers[pierIdx];
      auto pPier = bridge_desc.GetPier(pierIdx);
      if (0 < pierIdx)
      {
         auto spacing = GetPropertyList<Ifc4x3_add2, Ifc4x3_add2::IfcLengthMeasure>(pier, "pgsSpacing", "Back_Spacing");
         if (spacing.empty())
         {
            std::_tostringstream os;
            os << _T("Pier ") << LABEL_PIER(pierIdx) << _T(": Back_Spacing property in pgsSpacing property set not found");
            IFC_THROW(os.str().c_str());
         }

         auto girder_spacing = pPier->GetGirderSpacing(pgsTypes::Back);
         girder_spacing->SetMeasurementType(pgsTypes::MeasurementType::NormalToItem); // this is how the spacing is defined in the exporter
         girder_spacing->SetMeasurementLocation(pgsTypes::MeasurementLocation::AtCenterlineBearing);
         girder_spacing->ExpandAll();
         IndexType idx = 0;
         for (auto s : spacing)
         {
            girder_spacing->SetGirderSpacing(idx++, *s);
         }
      }

      if (pierIdx < nPiers - 1)
      {
         auto spacing = GetPropertyList<Ifc4x3_add2, Ifc4x3_add2::IfcLengthMeasure>(pier, "pgsSpacing", "Ahead_Spacing");
         if (spacing.empty())
         {
            std::_tostringstream os;
            os << _T("Pier ") << LABEL_PIER(pierIdx) << _T(": Ahead_Spacing property in pgsSpacing property set not found");
            IFC_THROW(os.str().c_str());
         }

         auto girder_spacing = pPier->GetGirderSpacing(pgsTypes::Ahead);
         girder_spacing->SetMeasurementType(pgsTypes::MeasurementType::NormalToItem); // this is how the spacing is defined in the exporter
         girder_spacing->SetMeasurementLocation(pgsTypes::MeasurementLocation::AtCenterlineBearing);

         girder_spacing->ExpandAll();
         IndexType idx = 0;
         for (auto s : spacing)
         {
            girder_spacing->SetGirderSpacing(idx++, *s);
         }
      }
   }

   SetGirderProperties(file, bridge_desc);

   ImportSlab(file, bridge_desc);

   pIBridgeDesc->SetBridgeDescription(bridge_desc);

   return CIfcImporter::ImportResult::Success;
}



void CIfcBridgeImporter::SetGirderProperties(IfcParse::IfcFile& file, CBridgeDescription2& bridge_desc)
{
   USES_CONVERSION;

   auto beams = file.instances_by_type<Ifc4x3_add2::IfcBeam>();
   Ifc4x3_add2::IfcBeam::list::ptr prestressed_beams(new Ifc4x3_add2::IfcBeam::list);
   for (auto beam : *beams)
   {
      auto predefined_type = GetPredefinedType<Ifc4x3_add2::IfcBeam, Ifc4x3_add2::IfcBeamType, Ifc4x3_add2::IfcBeamTypeEnum::Value>(beam);
      if (predefined_type.value_or(Ifc4x3_add2::IfcBeamTypeEnum::IfcBeamType_NOTDEFINED) == Ifc4x3_add2::IfcBeamTypeEnum::IfcBeamType_BEAM && HasClassification<Ifc4x3_add2>(beam, "usBridge_GirderPrestressedConcrete"))
         prestressed_beams->push(beam);
   }

   int beam_type_count = GetBeamTypeCount(file);

   bridge_desc.UseSameGirderForEntireBridge(beam_type_count == 1 ? true : false);
   GET_IFACE2(m_Importer.GetBroker(),ILibrary, pLibrary);
   if (bridge_desc.UseSameGirderForEntireBridge())
   {
      auto type = GetType<Ifc4x3_add2::IfcBeamType>(*beams->begin());
      auto girder_name = type->Name().get_value_or(std::string("Unknown"));
      auto girder_library_entry = pLibrary->GetGirderEntry(A2T(girder_name.c_str()));
      if (!girder_library_entry)
      {
         std::ostringstream os;
         os << "Girder type \"" << girder_name << "\" not found in the library";
         IFC_THROW(A2T(os.str().c_str()));
      }
      bridge_desc.SetGirderLibraryEntry(girder_library_entry);
      bridge_desc.SetGirderFamilyName(girder_library_entry->GetGirderFamilyName().c_str());

#pragma Reminder("WORKING HERE - This is assuming the first supported orientation. The IFC file doesn't have this information.")
      auto factory = girder_library_entry->GetBeamFactory();
      auto orientations = factory->GetSupportedGirderOrientation();
      bridge_desc.SetGirderOrientation(orientations.front());
   }

   for (auto beam : *prestressed_beams)
   {
      auto value = GetProperty<Ifc4x3_add2, Ifc4x3_add2::IfcLabel>(beam, "Pset_PrecastConcreteElementGeneral", "DesignLocationNumber");
      if (!value)
      {
         IFC_THROW(_T("DesignLocationNumber in Pset_PrecastConcreteElement property set not found"));
      }

      auto [spanIdx, gdrIdx] = ExtractSpanAndGirder(*value);
      auto fci = GetProperty<Ifc4x3_add2, Ifc4x3_add2::IfcPressureMeasure>(beam, "Pset_PrecastConcreteElementGeneral", "ReleaseStrength");
      if (!fci)
      {
         IFC_THROW(_T("ReleaseStrength property in Pset_PrecastConcreteElementGeneral property set not found"));
      }
      bridge_desc.GetGirderGroup(spanIdx)->GetGirder(gdrIdx)->GetSegment(0)->Material.Concrete.Fci = *fci;

      auto material = GetMaterial<Ifc4x3_add2>(beam);
      auto fc = GetMaterialProperty<Ifc4x3_add2, Ifc4x3_add2::IfcPressureMeasure>(material, "Pset_MaterialConcrete", "CompressiveStrength");
      if (!fc)
      {
         IFC_THROW(_T("CompressiveStrength property in Pset_PrecastConcreteElementGeneral property set not found"));
      }
      bridge_desc.GetGirderGroup(spanIdx)->GetGirder(gdrIdx)->GetSegment(0)->Material.Concrete.Fc = *fc;

      if (!bridge_desc.UseSameGirderForEntireBridge())
      {
         auto type = GetType<Ifc4x3_add2::IfcBeamType>(beam);
         auto girder_name = type->Name().get_value_or(std::string("Unknown"));
         auto girder_library_entry = pLibrary->GetGirderEntry(A2T(girder_name.c_str()));
         if (!girder_library_entry)
         {
            std::ostringstream os;
            os << "Girder type \"" << girder_name << "\" not found in the library";
            IFC_THROW(A2T(os.str().c_str()));
         }
         bridge_desc.GetGirderGroup(spanIdx)->GetGirder(gdrIdx)->SetGirderLibraryEntry(girder_library_entry);
         bridge_desc.SetGirderFamilyName(girder_library_entry->GetGirderFamilyName().c_str());

#pragma Reminder("WORKING HERE - This is assuming the first supported orientation. The IFC file doesn't have this information.")
#pragma Reminder("WORKING HERE - This is duplicate code from above. Simplify")
         auto factory = girder_library_entry->GetBeamFactory();
         auto orientations = factory->GetSupportedGirderOrientation();
         bridge_desc.SetGirderOrientation(orientations.front());
      }
   }
}

bool CIfcBridgeImporter::IsValidBridge(IfcParse::IfcFile& file, Ifc4x3_add2::IfcBridge* bridge)
{
   // must be a girder bridge
   if (bridge->PredefinedType().value_or(Ifc4x3_add2::IfcBridgeTypeEnum::IfcBridgeType_NOTDEFINED) != Ifc4x3_add2::IfcBridgeTypeEnum::IfcBridgeType_GIRDER)
      return false;

   // This check could be far less strict if we can assume the user provided us a PSG bridge.
   // If we can go from the girder line geometry and the girder name, mapped to a library entry,
   // that might be enough to actually do some work
   if (!HasValidGirders(file, bridge))
      return false;

   return true;
}

bool CIfcBridgeImporter::HasValidGirders(IfcParse::IfcFile& file, Ifc4x3_add2::IfcBridge* bridge)
{
   if (HasValidGirdersByTPF(file, bridge))
      return true;

   if (HasValidGirdersByOther(file, bridge))
      return true;

   m_Importer.AddNote(_T("One or more beams in the superstructure could not be identified as precast, prestressed concrete."));

   return false;
}

bool CIfcBridgeImporter::HasValidGirdersByTPF(IfcParse::IfcFile& file, Ifc4x3_add2::IfcBridge* bridge)
{
   auto superstructure = GetBridgePart(file, Ifc4x3_add2::IfcBridgePartTypeEnum::IfcBridgePartType_SUPERSTRUCTURE);

   auto beams = file.instances_by_type<Ifc4x3_add2::IfcBeam>();
   bool valid_beams = true;
   for (auto beam : *beams)
   {
      // beam must be contained in the spatial structure of the superstructure
      auto related_elements = beam->ContainedInStructure();
      if (related_elements)
      {
         for (auto related_element : *related_elements)
         {
            if (related_element->RelatingStructure() != superstructure)
               continue;
         }
      }
      else
      {
         continue; // not in a spatial structure
      }

      // beam must be IfcBeam.BEAM
      auto predefined_type = GetPredefinedType<Ifc4x3_add2::IfcBeam, Ifc4x3_add2::IfcBeamType, Ifc4x3_add2::IfcBeamTypeEnum::Value>(beam);
      if (predefined_type.value_or(Ifc4x3_add2::IfcBeamTypeEnum::IfcBeamType_NOTDEFINED) == Ifc4x3_add2::IfcBeamTypeEnum::IfcBeamType_BEAM)
      {
         // beams must be classified as precast girders
         auto assembly_place = GetPropertyEnum<Ifc4x3_add2, Ifc4x3_add2::IfcLabel>(beam, "Pset_ConcreteElementGeneral", "AssemblyPlace");
         auto casting_method = GetPropertyEnum<Ifc4x3_add2, Ifc4x3_add2::IfcLabel>(beam, "Pset_ConcreteElementGeneral", "CastingMethod");

         // level 1 check - if there is an assembly_place and casting_method, use that for the determination
         // otherwise use the more specific usBridge classification
         if (assembly_place && casting_method)
         {
            bool is_factory_assembled = (std::string(*assembly_place) == std::string("FACTORY"));
            bool is_precast = (std::string(*casting_method) == std::string("PRECAST"));

            if (!is_factory_assembled || !is_precast)
            {
               valid_beams = false;
               break;
            }
         }
         else
         {
            // This is not a great requirement - TPF has decided that this is the only way to determine if a beam is prestressed concrete
            // I think we are going to find that many don't use this classification
            if (!HasClassification<Ifc4x3_add2>(beam, "usBridge_GirderPrestressedConcrete"))
            {
               valid_beams = false;
               break;
            }
         }
      }
   }

   return valid_beams;
}

bool CIfcBridgeImporter::HasValidGirdersByOther(IfcParse::IfcFile& file, Ifc4x3_add2::IfcBridge* bridge)
{
   // This function attempts to determine if all the beams in the superstructure are precast concrete.
   // The basic idea is that the beams are IfcElementAssembly and they are factory assembled girders.
   // A factory assembled girder alone is not enough to claim the girders are precast concrete.
   //
   // This function may be a bad idea... keep it for now, but be skeptical
   auto superstructure = GetBridgePart(file, Ifc4x3_add2::IfcBridgePartTypeEnum::IfcBridgePartType_SUPERSTRUCTURE);

   auto element_assemblies = file.instances_by_type<Ifc4x3_add2::IfcElementAssembly>();
   bool valid_beams = true;
   for (auto element_assembly : *element_assemblies)
   {
      // beam must be contained in the spatial structure of the superstructure
      auto related_elements = element_assembly->ContainedInStructure();
      if (related_elements)
      {
         for (auto related_element : *related_elements)
         {
            if (related_element->RelatingStructure() != superstructure)
               continue;
         }
      }
      else
      {
         continue; // not in a spatial structure
      }

      auto assembly_place = element_assembly->AssemblyPlace();
      bool is_factory_assembled = (assembly_place.value_or(Ifc4x3_add2::IfcAssemblyPlaceEnum::IfcAssemblyPlace_NOTDEFINED) == Ifc4x3_add2::IfcAssemblyPlaceEnum::IfcAssemblyPlace_FACTORY);

      auto predefined_type = GetPredefinedType<Ifc4x3_add2::IfcElementAssembly, Ifc4x3_add2::IfcElementAssemblyType, Ifc4x3_add2::IfcElementAssemblyTypeEnum::Value>(element_assembly);
      bool is_girder = (predefined_type.value_or(Ifc4x3_add2::IfcElementAssemblyTypeEnum::IfcElementAssemblyType_NOTDEFINED) == Ifc4x3_add2::IfcElementAssemblyTypeEnum::IfcElementAssemblyType_GIRDER);
      if (!is_factory_assembled || !is_girder)
      {
         valid_beams = false;
         break;
      }
   }

   return valid_beams;
}

Ifc4x3_add2::IfcBridge* CIfcBridgeImporter::GetBridge(IfcParse::IfcFile& file)
{
   USES_CONVERSION;

   auto bridges = file.instances_by_type<Ifc4x3_add2::IfcBridge>();

   if (1 <= bridges->size())
   {
      std::vector<Ifc4x3_add2::IfcBridge*> valid_bridges;
      for (auto bridge : *bridges)
      {
         if (IsValidBridge(file, bridge))
            valid_bridges.push_back(bridge);
      }

      if (valid_bridges.size() == 0)
      {
         m_Importer.AddNote(_T("IFC model does not contain a bridge that is compatible with this software."));
      }
      else
      {
         int result = 0;
         if (1 < valid_bridges.size())
         {
            // only prompt if there is more than one bridge
            std::ostringstream os;
            for (auto bridge : valid_bridges)
            {
               auto strLabel = (bridge->Name() ? *(bridge->Name()) : bridge->Description() ? *(bridge->Description()) : "Unnamed");
               os << strLabel << std::endl;
            }

            auto result = AfxChoose(_T("Select Bridge"), _T("Select bridge to import"), A2T(os.str().c_str()), 0, TRUE);
            if (result < 0)
               return nullptr;
         }

         return valid_bridges[result];
      }
   }

   return nullptr;
}

void CIfcBridgeImporter::ImportSlab(IfcParse::IfcFile& file, CBridgeDescription2& bridge_desc)
{
   auto slabs = file.instances_by_type<Ifc4x3_add2::IfcSlab>();
   auto it = std::find_if(slabs->begin(), slabs->end(), [](const auto& slab) {return slab->PredefinedType() == Ifc4x3_add2::IfcSlabTypeEnum::IfcSlabType_FLOOR; });
   if (it == slabs->end())
      return; // no slabs

   auto slab = *it;

   auto stations = GetPropertyList<Ifc4x3_add2, Ifc4x3_add2::IfcLengthMeasure>(slab, "pgsDeck", "Stations");
   if (stations.empty())
   {
      IFC_THROW(_T("Stations property in pgsDeck property set not found"));
   }

   auto left_edges = GetPropertyList<Ifc4x3_add2, Ifc4x3_add2::IfcLengthMeasure>(slab, "pgsDeck", "LeftEdges");
   if (left_edges.empty())
   {
      IFC_THROW(_T("LeftEdges property in pgsDeck property set not found"));
   }

   auto right_edges = GetPropertyList<Ifc4x3_add2, Ifc4x3_add2::IfcLengthMeasure>(slab, "pgsDeck", "RightEdges");
   if (right_edges.empty())
   {
      IFC_THROW(_T("RightEdges property in pgsDeck property set not found"));
   }

   if (stations.size() != left_edges.size() || stations.size() != right_edges.size())
   {
      IFC_THROW(_T("Stations, LeftEdges, and RightEdges properties in pgsDeck property set must have the same number of values"));
   }

   auto* pDeck = bridge_desc.GetDeckDescription();

   auto* gross_depth = GetProperty<Ifc4x3_add2, Ifc4x3_add2::IfcLengthMeasure>(slab, "pgsDeck", "GrossDepth");
   if (gross_depth)
   {
      pDeck->GrossDepth = *gross_depth;
   }
   else
   {
      IFC_THROW(_T("GrossDepth property in pgsDeck property set not found"));
   }

   auto* left_edge_depth = GetProperty<Ifc4x3_add2, Ifc4x3_add2::IfcLengthMeasure>(slab, "pgsDeck", "LeftEdgeDepth");
   if (left_edge_depth)
   {
      pDeck->OverhangEdgeDepth[pgsTypes::stLeft] = *left_edge_depth;
   }
   else
   {
      IFC_THROW(_T("LeftEdgeDepth property in pgsDeck property set not found"));
   }

   auto* right_edge_depth = GetProperty<Ifc4x3_add2, Ifc4x3_add2::IfcLengthMeasure>(slab, "pgsDeck", "RightEdgeDepth");
   if (right_edge_depth)
   {
      pDeck->OverhangEdgeDepth[pgsTypes::stRight] = *right_edge_depth;
   }
   else
   {
      IFC_THROW(_T("RightEdgeDepth property in pgsDeck property set not found"));
   }

   pDeck->DeckEdgePoints.clear();
   for (auto&& [station, left_edge, right_edge] : boost::combine(stations, left_edges, right_edges))
   {
      CDeckPoint deck_point;

      // dummy, default values
      deck_point.MeasurementType = pgsTypes::OffsetMeasurementType::omtBridge;
      deck_point.LeftTransitionType = stations.size() == 1 ? pgsTypes::DeckPointTransitionType::dptParallel : pgsTypes::DeckPointTransitionType::dptLinear;
      deck_point.RightTransitionType = stations.size() == 1 ? pgsTypes::DeckPointTransitionType::dptParallel : pgsTypes::DeckPointTransitionType::dptLinear;

      deck_point.Station = *station;
      deck_point.LeftEdge = *left_edge;
      deck_point.RightEdge = *right_edge;
      pDeck->DeckEdgePoints.emplace_back(deck_point);
   }
}
