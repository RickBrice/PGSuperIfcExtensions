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
#include "IfcBridgeImporter.h"
#include "IfcImporter.h"
#include "Properties.h"
#include "IfcImporterException.h"
#include "USBridge_Classifications.h"
#include "Utilities.h"

#include "BeamSpacing.h"
#include "DeckSlab.h"
#include "Piers.h"

#include <MFCTools\Prompts.h>
#include <boost/range/combine.hpp>
#include <psgLib/BridgeDescription2.h>
#include <psgLib/GirderLabel.h>
#include <psgLib/GirderLibraryEntry.h>



CIfcBridgeImporter::CIfcBridgeImporter(CIfcImporter& importer) :
   m_Importer(importer)
{
}

CIfcImporter::ImportResult CIfcBridgeImporter::Import(IfcParse::IfcFile& file, bool bDeriveAlignmentFromDeck)
{
   auto bridge = GetBridge(file);
   if (bridge == nullptr)
      return CIfcImporter::ImportResult::NotFound;

   if (bDeriveAlignmentFromDeck)
   {
      if (!DeriveAlignmentFromDeck(file))
      {
         WBFL::System::Logger::Info(_T("Failed to derive alignment from deck slab."));
         //IFC_THROW(_T("Failed to derive alignment from deck slab."));
         return CIfcImporter::ImportResult::Fail;
      }
   }

   auto nPiers = get_pier_count(file);

   SpanIndexType nSpans = INVALID_INDEX;
   auto value = GetProperty<IfcSchema, IfcSchema::IfcInteger>(bridge, "usBridge_BridgeCommon", "usBridge_NumberOfSpans");
   if (value)
   {
      nSpans = (SpanIndexType)(*value);
      if (nSpans != nPiers - 1)
         IFC_THROW(_T("Number of spans modeled does not match number of spans in usBridge_BridgeCommon property set"));
   }
   else
   {
      WBFL::System::Logger::Info(_T("usBridge_NumberOfSpans property not found in usBridge_BridgeCommon property set"));
      nSpans = nPiers - 1; // derive number of spans from nPiers so we don't rely on custom property set
      //IFC_THROW(_T("usBridge_NumberOfSpans property not found in usBridge_BridgeCommon property set"));
   }

   std::vector<GirderIndexType> nGirders;
   nGirders.assign(nSpans, 0);
   auto beam_ids = get_beam_ids(file);
   IndexType girder_count = beam_ids.size();
   IndexType girders_per_span = girder_count / nSpans;
   IndexType girders_processed = 0;
   for (auto id : beam_ids)
   {
      auto beam = file.instance_by_id(id)->as<IfcSchema::IfcBeam>();
      auto girder_key = get_girder_key(beam);
      if (girder_key == CGirderKey())
      {
         WBFL::System::Logger::Info(_T("Using assumed girder key."));
         girder_key.groupIndex = girders_processed / girders_per_span; // assume girders are evenly distributed across spans
         //IFC_THROW(_T("DesignLocationNumber property not found in Pset_PrecastConcreteElementGeneral"));
      }
      nGirders[girder_key.groupIndex]++;
      girders_processed++;
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
   }


   if (bSameNumGirdersInAllSpans)
   {
      bridge_desc.UseSameNumberOfGirdersInAllGroups(true);
      bridge_desc.SetGirderCount(nGirders.front());
   }
   else 
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

   // Per TPF modeling guidance, piers and abutments are different types.
   // Get the abutments and piers and put into a single vector because we need to treat them the same in PGSuper
   std::vector<IfcSchema::IfcBridgePart*> abutments = GetBridgeParts(file, IfcSchema::IfcBridgePartTypeEnum::Value::IfcBridgePartType_ABUTMENT);
   std::vector<IfcSchema::IfcBridgePart*> piers = GetBridgeParts(file, IfcSchema::IfcBridgePartTypeEnum::Value::IfcBridgePartType_PIER);
   piers.insert(piers.begin(), abutments.front());
   piers.insert(piers.end(), abutments.back());

   ASSERT(nPiers == bridge_desc.GetPierCount());
   std::set<double> stations; // we don't know the order the piers are defined in the model... could be anything. Set will sort stations
   for (PierIndexType pierIdx = 0; pierIdx < nPiers; pierIdx++)
   {
      // get the pier station from the positioning element and set it on the PGSuper pier
      auto pier = piers[pierIdx];
      double station = get_pier_station(m_Importer.GetBroker(),file, pierIdx, pier);
      stations.insert(station);
   }

   PierIndexType pierIdx = 0;
   for (auto& station : stations)
   {
      auto pPier = bridge_desc.GetPier(pierIdx++);
      pPier->SetStation(station);
   }

   // This code is commented out, because the property set is no longer used. Spacing
   // is computed from the bridge model geometry. See below
   // 
   //// NOTE: Girder spacing is stored in a property set attached to IfcPier. It is in its own
   //// property set pgsSpacing so it doesn't interfere with the stationing property.
   //bridge_desc.SetGirderSpacingType(pgsTypes::SupportedBeamSpacing::sbsGeneral);
   //for (PierIndexType pierIdx = 0; pierIdx < nPiers; pierIdx++)
   //{
   //   auto pier = piers[pierIdx];
   //   auto pPier = bridge_desc.GetPier(pierIdx);
   //   if (0 < pierIdx)
   //   {
   //      auto spacing = GetPropertyList<IfcSchema, IfcSchema::IfcLengthMeasure>(pier, "pgsSpacing", "Back_Spacing");
   //      if (spacing.empty())
   //      {
   //         std::_tostringstream os;
   //         os << _T("Pier ") << LABEL_PIER(pierIdx) << _T(": Back_Spacing property in pgsSpacing property set not found");
   //         IFC_THROW(os.str().c_str());
   //      }

   //      auto girder_spacing = pPier->GetGirderSpacing(pgsTypes::Back);
   //      girder_spacing->SetMeasurementType(pgsTypes::MeasurementType::NormalToItem); // this is how the spacing is defined in the exporter
   //      girder_spacing->SetMeasurementLocation(pgsTypes::MeasurementLocation::AtCenterlineBearing);
   //      girder_spacing->ExpandAll();
   //      IndexType idx = 0;
   //      for (auto s : spacing)
   //      {
   //         girder_spacing->SetGirderSpacing(idx++, *s);
   //      }
   //   }

   //   if (pierIdx < nPiers - 1)
   //   {
   //      auto spacing = GetPropertyList<IfcSchema, IfcSchema::IfcLengthMeasure>(pier, "pgsSpacing", "Ahead_Spacing");
   //      if (spacing.empty())
   //      {
   //         std::_tostringstream os;
   //         os << _T("Pier ") << LABEL_PIER(pierIdx) << _T(": Ahead_Spacing property in pgsSpacing property set not found");
   //         IFC_THROW(os.str().c_str());
   //      }

   //      auto girder_spacing = pPier->GetGirderSpacing(pgsTypes::Ahead);
   //      girder_spacing->SetMeasurementType(pgsTypes::MeasurementType::NormalToItem); // this is how the spacing is defined in the exporter
   //      girder_spacing->SetMeasurementLocation(pgsTypes::MeasurementLocation::AtCenterlineBearing);

   //      girder_spacing->ExpandAll();
   //      IndexType idx = 0;
   //      for (auto s : spacing)
   //      {
   //         girder_spacing->SetGirderSpacing(idx++, *s);
   //      }
   //   }
   //}

   // obtain beam spacing from the bridge model geometry
   auto [start_spacing, end_spacing] = get_beam_spacing(m_Importer.GetBroker(),file);

   // assume general spacing for now, but in the future analyze the spacing
   // data and see if it is the same for the entire bridge, same for a span, or girder by girder
   bridge_desc.SetGirderSpacingType(pgsTypes::SupportedBeamSpacing::sbsGeneral);
   for( auto spanIdx = 0; spanIdx < nSpans; spanIdx++)
   {
      auto pSpan = bridge_desc.GetSpan(spanIdx);
      auto start_girder_spacing = pSpan->GetPrevPier()->GetGirderSpacing(pgsTypes::Ahead);
      start_girder_spacing->SetMeasurementType(pgsTypes::MeasurementType::AlongItem);
      // Right now the beam spacing is based on the girder end points, so AtPierLine is not exactly correct.
      // As work with extracting key parameters from the model geometry progresses, eventually the pier
      // line geometry will be found, and then the girder spacing will be measured at the pier line or CL Bearing.
      // See also below for end of girder spacing.
      start_girder_spacing->SetMeasurementLocation(pgsTypes::MeasurementLocation::AtPierLine);
      start_girder_spacing->ExpandAll();
      IndexType idx = 0;
      for (auto s : start_spacing[spanIdx])
      {
         start_girder_spacing->SetGirderSpacing(idx++, s);
      }

      auto end_girder_spacing = pSpan->GetNextPier()->GetGirderSpacing(pgsTypes::Back);
      end_girder_spacing->SetMeasurementType(pgsTypes::MeasurementType::AlongItem);
      end_girder_spacing->SetMeasurementLocation(pgsTypes::MeasurementLocation::AtPierLine);
      end_girder_spacing->ExpandAll();
      idx = 0;
      for (auto s : end_spacing[spanIdx])
      {
         end_girder_spacing->SetGirderSpacing(idx++, s);
      }
   }

   SetGirderProperties(file, bridge_desc);

   ImportSlab(file, bridge_desc);

   pIBridgeDesc->SetBridgeDescription(bridge_desc);

   Experiment(file);

   return CIfcImporter::ImportResult::Success;
}



void CIfcBridgeImporter::SetGirderProperties(IfcParse::IfcFile& file, CBridgeDescription2& bridge_desc)
{
   USES_CONVERSION;

   auto beams = file.instances_by_type<IfcSchema::IfcBeam>();
   IfcSchema::IfcBeam::list::ptr prestressed_beams(new IfcSchema::IfcBeam::list);
   for (auto beam : *beams)
   {
      auto predefined_type = GetPredefinedType<IfcSchema::IfcBeam, IfcSchema::IfcBeamType, IfcSchema::IfcBeamTypeEnum::Value>(beam);
      if (predefined_type.value_or(IfcSchema::IfcBeamTypeEnum::IfcBeamType_NOTDEFINED) == IfcSchema::IfcBeamTypeEnum::IfcBeamType_BEAM && HasClassification<IfcSchema>(beam, "usBridge_GirderPrestressedConcrete"))
         prestressed_beams->push(beam);
   }

   if (prestressed_beams->size() == 0)
   {
      WBFL::System::Logger::Info("Did not find IfcBeam in the superstructure spatial structure classified as usBridge_GirderPrestressedConcrete. Assuming all superstructure IfcBeam.BEAM are prestressed girders.");
      auto beam_ids = get_beam_ids(file);
      for (auto id : beam_ids)
      {
         prestressed_beams->push(file.instance_by_id(id)->as<IfcSchema::IfcBeam>());
      }
   }

   int beam_type_count = GetBeamTypeCount(file);

   bridge_desc.UseSameGirderForEntireBridge(beam_type_count == 1 ? true : false);
   if (bridge_desc.UseSameGirderForEntireBridge())
   {
      auto girder_library_entry = GetGirderLibraryEntry(*beams->begin());
      bridge_desc.SetGirderLibraryEntry(girder_library_entry);
      bridge_desc.SetGirderFamilyName(girder_library_entry->GetGirderFamilyName().c_str());

#pragma Reminder("WORKING HERE - This is assuming the first supported orientation. The IFC file doesn't have this information.")
      // should get the orientation from the girder geometry and then try to match it with the supported orientations
      // if no match, then log warning and use the default.
      auto factory = girder_library_entry->GetBeamFactory();
      auto orientations = factory->GetSupportedGirderOrientation();
      bridge_desc.SetGirderOrientation(orientations.front());
   }

   IndexType girders_processed = 0;
   for (auto beam : *prestressed_beams)
   {
      auto girder_key = get_girder_key(beam);
      if (girder_key == CGirderKey())
      {
         WBFL::System::Logger::Info(_T("Using assumed girder key."));
         girder_key.groupIndex = 0;
         girder_key.girderIndex = girders_processed++;
      }

      auto fci = GetProperty<IfcSchema, IfcSchema::IfcPressureMeasure>(beam, "Pset_PrecastConcreteElementGeneral", "ReleaseStrength");
      if (fci)
      {
         bridge_desc.GetGirderGroup(girder_key.groupIndex)->GetGirder(girder_key.girderIndex)->GetSegment(0)->Material.Concrete.Fci = *fci;
      }
      else
      {
         WBFL::System::Logger::Info(_T("ReleaseStrength property in Pset_PrecastConcreteElementGeneral property set not found"));
         //IFC_THROW(_T("ReleaseStrength property in Pset_PrecastConcreteElementGeneral property set not found"));
      }

      auto material = GetMaterial<IfcSchema>(beam);
      if (material)
      {
         auto fc = GetMaterialProperty<IfcSchema, IfcSchema::IfcPressureMeasure>(material, "Pset_MaterialConcrete", "CompressiveStrength");
         if (fc)
         {
            bridge_desc.GetGirderGroup(girder_key.groupIndex)->GetGirder(girder_key.girderIndex)->GetSegment(0)->Material.Concrete.Fc = *fc;
         }
         else
         {
            WBFL::System::Logger::Info(_T("CompressiveStrength property in Pset_PrecastConcreteElementGeneral property set not found"));
            //IFC_THROW(_T("CompressiveStrength property in Pset_PrecastConcreteElementGeneral property set not found"));
         }
      }
      else
      {
         WBFL::System::Logger::Info(_T("Materials are not associated with the beam"));
         //IFC_THROW(_T("Materials are not associated with the beam"));
      }

      if (!bridge_desc.UseSameGirderForEntireBridge())
      {
         auto girder_library_entry = GetGirderLibraryEntry(beam);
         bridge_desc.GetGirderGroup(girder_key.groupIndex)->GetGirder(girder_key.girderIndex)->SetGirderLibraryEntry(girder_library_entry);
         bridge_desc.SetGirderFamilyName(girder_library_entry->GetGirderFamilyName().c_str());

#pragma Reminder("WORKING HERE - This is assuming the first supported orientation. The IFC file doesn't have this information.")
#pragma Reminder("WORKING HERE - This is duplicate code from above. Simplify")
         auto factory = girder_library_entry->GetBeamFactory();
         auto orientations = factory->GetSupportedGirderOrientation();
         bridge_desc.SetGirderOrientation(orientations.front());
      }
   }
}

bool CIfcBridgeImporter::IsValidBridge(IfcParse::IfcFile& file, IfcSchema::IfcBridge* bridge)
{
   // must be a girder bridge
   if (bridge->PredefinedType().value_or(IfcSchema::IfcBridgeTypeEnum::IfcBridgeType_NOTDEFINED) != IfcSchema::IfcBridgeTypeEnum::IfcBridgeType_GIRDER)
      return false;

   // This check could be far less strict if we can assume the user provided us a PSG bridge.
   // If we can go from the girder line geometry and the girder name, mapped to a library entry,
   // that might be enough to actually do some work
   //if (!HasValidGirders(file, bridge))
   //   return false;

   return true;
}

bool CIfcBridgeImporter::HasValidGirders(IfcParse::IfcFile& file, IfcSchema::IfcBridge* bridge)
{
   if (HasValidGirdersByTPF(file, bridge))
      return true;

   if (HasValidGirdersByOther(file, bridge))
      return true;

   WBFL::System::Logger::Info(_T("One or more beams in the superstructure could not be identified as precast, prestressed concrete."));

   return false;
}

bool CIfcBridgeImporter::HasValidGirdersByTPF(IfcParse::IfcFile& file, IfcSchema::IfcBridge* bridge)
{
   auto superstructure = GetBridgePart(file, IfcSchema::IfcBridgePartTypeEnum::IfcBridgePartType_SUPERSTRUCTURE);

   auto beams = file.instances_by_type<IfcSchema::IfcBeam>();
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
      auto predefined_type = GetPredefinedType<IfcSchema::IfcBeam, IfcSchema::IfcBeamType, IfcSchema::IfcBeamTypeEnum::Value>(beam);
      if (predefined_type.value_or(IfcSchema::IfcBeamTypeEnum::IfcBeamType_NOTDEFINED) == IfcSchema::IfcBeamTypeEnum::IfcBeamType_BEAM)
      {
         // beams must be classified as precast girders
         auto assembly_place = GetPropertyEnum<IfcSchema, IfcSchema::IfcLabel>(beam, "Pset_ConcreteElementGeneral", "AssemblyPlace");
         auto casting_method = GetPropertyEnum<IfcSchema, IfcSchema::IfcLabel>(beam, "Pset_ConcreteElementGeneral", "CastingMethod");

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
            if (!HasClassification<IfcSchema>(beam, "usBridge_GirderPrestressedConcrete"))
            {
               valid_beams = false;
               break;
            }
         }
      }
   }

   return valid_beams;
}

bool CIfcBridgeImporter::HasValidGirdersByOther(IfcParse::IfcFile& file, IfcSchema::IfcBridge* bridge)
{
   // This function attempts to determine if all the beams in the superstructure are precast concrete.
   // The basic idea is that the beams are IfcElementAssembly and they are factory assembled girders.
   // A factory assembled girder alone is not enough to claim the girders are precast concrete.
   //
   // This function may be a bad idea... keep it for now, but be skeptical
   auto superstructure = GetBridgePart(file, IfcSchema::IfcBridgePartTypeEnum::IfcBridgePartType_SUPERSTRUCTURE);

   auto element_assemblies = file.instances_by_type<IfcSchema::IfcElementAssembly>();
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
      bool is_factory_assembled = (assembly_place.value_or(IfcSchema::IfcAssemblyPlaceEnum::IfcAssemblyPlace_NOTDEFINED) == IfcSchema::IfcAssemblyPlaceEnum::IfcAssemblyPlace_FACTORY);

      auto predefined_type = GetPredefinedType<IfcSchema::IfcElementAssembly, IfcSchema::IfcElementAssemblyType, IfcSchema::IfcElementAssemblyTypeEnum::Value>(element_assembly);
      bool is_girder = (predefined_type.value_or(IfcSchema::IfcElementAssemblyTypeEnum::IfcElementAssemblyType_NOTDEFINED) == IfcSchema::IfcElementAssemblyTypeEnum::IfcElementAssemblyType_GIRDER);
      if (!is_factory_assembled || !is_girder)
      {
         valid_beams = false;
         break;
      }
   }

   return valid_beams;
}

IfcSchema::IfcBridge* CIfcBridgeImporter::GetBridge(IfcParse::IfcFile& file)
{
   USES_CONVERSION;

   auto bridges = file.instances_by_type<IfcSchema::IfcBridge>();

   if (1 <= bridges->size())
   {
      std::vector<IfcSchema::IfcBridge*> valid_bridges;
      for (auto bridge : *bridges)
      {
         if (IsValidBridge(file, bridge))
            valid_bridges.push_back(bridge);
      }

      if (valid_bridges.size() == 0)
      {
         WBFL::System::Logger::Info(_T("IFC model does not contain a bridge that is compatible with this software."));
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
   auto slabs = file.instances_by_type<IfcSchema::IfcSlab>();
   auto it = std::find_if(slabs->begin(), slabs->end(), [](const auto& slab) {return slab->PredefinedType() == IfcSchema::IfcSlabTypeEnum::IfcSlabType_FLOOR; });
   if (it == slabs->end())
      return; // no slabs

   auto slab = *it;

   auto deck_points = get_deck_slab(m_Importer.GetBroker(), file);

   //auto stations = GetPropertyList<IfcSchema, IfcSchema::IfcLengthMeasure>(slab, "pgsDeck", "Stations");
   //if (stations.empty())
   //{
   //   IFC_THROW(_T("Stations property in pgsDeck property set not found"));
   //}

   //auto left_edges = GetPropertyList<IfcSchema, IfcSchema::IfcLengthMeasure>(slab, "pgsDeck", "LeftEdges");
   //if (left_edges.empty())
   //{
   //   IFC_THROW(_T("LeftEdges property in pgsDeck property set not found"));
   //}

   //auto right_edges = GetPropertyList<IfcSchema, IfcSchema::IfcLengthMeasure>(slab, "pgsDeck", "RightEdges");
   //if (right_edges.empty())
   //{
   //   IFC_THROW(_T("RightEdges property in pgsDeck property set not found"));
   //}

   //if (stations.size() != left_edges.size() || stations.size() != right_edges.size())
   //{
   //   IFC_THROW(_T("Stations, LeftEdges, and RightEdges properties in pgsDeck property set must have the same number of values"));
   //}

   auto* pDeck = bridge_desc.GetDeckDescription();

   // Use default for now, but ultimately need to cut a cross section through the deck to get these parameters
   //auto* gross_depth = GetProperty<IfcSchema, IfcSchema::IfcLengthMeasure>(slab, "pgsDeck", "GrossDepth");
   //if (gross_depth)
   //{
   //   pDeck->GrossDepth = *gross_depth;
   //}
   //else
   //{
   //   IFC_THROW(_T("GrossDepth property in pgsDeck property set not found"));
   //}

   //auto* left_edge_depth = GetProperty<IfcSchema, IfcSchema::IfcLengthMeasure>(slab, "pgsDeck", "LeftEdgeDepth");
   //if (left_edge_depth)
   //{
   //   pDeck->OverhangEdgeDepth[pgsTypes::stLeft] = *left_edge_depth;
   //}
   //else
   //{
   //   IFC_THROW(_T("LeftEdgeDepth property in pgsDeck property set not found"));
   //}

   //auto* right_edge_depth = GetProperty<IfcSchema, IfcSchema::IfcLengthMeasure>(slab, "pgsDeck", "RightEdgeDepth");
   //if (right_edge_depth)
   //{
   //   pDeck->OverhangEdgeDepth[pgsTypes::stRight] = *right_edge_depth;
   //}
   //else
   //{
   //   IFC_THROW(_T("RightEdgeDepth property in pgsDeck property set not found"));
   //}

   // building deck based on data in custom psets
   //pDeck->DeckEdgePoints.clear();
   //for (auto&& [station, left_edge, right_edge] : boost::combine(stations, left_edges, right_edges))
   //{
   //   CDeckPoint deck_point;

   //   // dummy, default values
   //   deck_point.MeasurementType = pgsTypes::OffsetMeasurementType::omtBridge;
   //   deck_point.LeftTransitionType = stations.size() == 1 ? pgsTypes::DeckPointTransitionType::dptParallel : pgsTypes::DeckPointTransitionType::dptLinear;
   //   deck_point.RightTransitionType = stations.size() == 1 ? pgsTypes::DeckPointTransitionType::dptParallel : pgsTypes::DeckPointTransitionType::dptLinear;

   //   deck_point.Station = *station;
   //   deck_point.LeftEdge = *left_edge;
   //   deck_point.RightEdge = *right_edge;
   //   pDeck->DeckEdgePoints.emplace_back(deck_point);
   //}

   pDeck->DeckEdgePoints.clear();
   for (auto&& [station, edge] : deck_points)
   {
      auto [left_edge, right_edge] = edge;
      CDeckPoint deck_point;

      // dummy, default values
      deck_point.MeasurementType = pgsTypes::OffsetMeasurementType::omtAlignment; // deck geometry from IFC file is computed relative to the alignment
      deck_point.LeftTransitionType = deck_points.size() == 1 ? pgsTypes::DeckPointTransitionType::dptParallel : pgsTypes::DeckPointTransitionType::dptLinear;
      deck_point.RightTransitionType = deck_points.size() == 1 ? pgsTypes::DeckPointTransitionType::dptParallel : pgsTypes::DeckPointTransitionType::dptLinear;

      deck_point.Station = station;
      deck_point.LeftEdge = fabs(left_edge);
      deck_point.RightEdge = fabs(right_edge);
      pDeck->DeckEdgePoints.emplace_back(deck_point);
   }
}

const GirderLibraryEntry* CIfcBridgeImporter::GetGirderLibraryEntry(IfcSchema::IfcBeam* beam)
{
   USES_CONVERSION;

   GET_IFACE2(m_Importer.GetBroker(), ILibrary, pLibrary);

   auto type = GetType<IfcSchema::IfcBeamType>(beam);
   auto girder_name = type ? type->Name().get_value_or(std::string("Unknown")) : std::string("Unknown");
   auto girder_library_entry = pLibrary->GetGirderEntry(A2T(girder_name.c_str()));
   if (!girder_library_entry)
   {
      // Matching girder type in the library is a bad implementation.
      // Should be creating a new library entry for this girder type.
      // Could not match the girder type. So the program doesn't crap out,
      // get the first I-Beam type and substitute it. Not a great solution,
      // but the current focus is loading files that don't conform to the
      // AASHTO/TPF data standards. We want to be able to load any model
      // with a precast girder beam bridge
      GET_IFACE2(m_Importer.GetBroker(), ILibraryNames, pLibNames);
      std::vector<std::_tstring> names;
      pLibNames->EnumGirderNames(_T("I-Beam"), &names); // huge assumption that we are dealing with I beams.
      auto substitute_girder_name = names.front();

      girder_library_entry = pLibrary->GetGirderEntry(substitute_girder_name.c_str());

      std::ostringstream os;
      beam->toString(os);
      os << std::endl;
      os << "Girder type \"" << girder_name << "\" not found in the library, substituting " << T2A(substitute_girder_name.c_str());

      WBFL::System::Logger::Info(os.str().c_str());

      WBFL::System::Logger::Info(A2T(os.str().c_str()));
   }
   return girder_library_entry;
}

bool CIfcBridgeImporter::DeriveAlignmentFromDeck(IfcParse::IfcFile& file)
{
   return create_alignment_from_deck(m_Importer.GetBroker(), file);
}

void CIfcBridgeImporter::Experiment(IfcParse::IfcFile& file)
{
   //get_beam_spacing(file);
}

