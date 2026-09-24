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

#include <numeric>

#include <MFCTools\Prompts.h>
#include <boost/range/combine.hpp>
#include <psgLib/BridgeDescription2.h>
#include <psgLib/GirderLabel.h>
#include <psgLib/GirderLibraryEntry.h>



CIfcBridgeImporter::CIfcBridgeImporter(CIfcImporter& importer) :
   m_Importer(importer)
{
}

CIfcImporter::ImportResult CIfcBridgeImporter::Import(ifcopenshell::file& file, bool bDeriveAlignmentFromDeck)
{
   WBFL::System::Logger::Info(_T("Importing from bridge IFC file."));
   auto bridge = GetBridge(file);
   if (!bridge)
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
   auto value = GetProperty<IfcSchema, IfcSchema::IfcInteger>(bridge, "usBrPset_BridgeGeometry", "NumberOfSpans");
   if (value)
   {
      nSpans = (SpanIndexType)(int64_t)(*value);
      if (nSpans != nPiers - 1)
         IFC_THROW(_T("Number of spans modeled does not match number of spans in usBrPset_BridgeGeometry property set"));
   }
   else
   {
      WBFL::System::Logger::Info(_T("NumberOfSpans property not found in usBrPset_BridgeGeometry property set"));
      nSpans = nPiers - 1; // derive number of spans from nPiers so we don't rely on custom property set
      //IFC_THROW(_T("usBridge_NumberOfSpans property not found in usBrPset_BridgeGeometry property set"));
   }

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
   }

   //
   // Position the abutments and piers
   //

   // Per TPF modeling guidance, piers and abutments are different types.
   // Get the abutments and piers and put into a single vector because we need to treat them the same in PGSuper.
   // The order doesn't matter, the piers are sorted by station
   std::vector<IfcSchema::IfcBridgePart> piers = GetBridgeParts(file, IfcSchema::IfcBridgePartTypeEnum::Value::IfcBridgePartType_ABUTMENT);
   auto interior_piers = GetBridgeParts(file, IfcSchema::IfcBridgePartTypeEnum::Value::IfcBridgePartType_PIER);
   piers.insert(piers.end(), interior_piers.begin(), interior_piers.end());

   ASSERT(nPiers == piers.size());
   ASSERT(nPiers == bridge_desc.GetPierCount());
   // station and direction (if the model has it) of each pier, sorted by station because we don't know the order the piers are defined in the model
   std::vector<std::pair<Float64, std::optional<Eigen::Vector2d>>> pier_locations;
   for (PierIndexType pierIdx = 0; pierIdx < nPiers; pierIdx++)
   {
      // get the pier station from the positioning element
      pier_locations.emplace_back(get_pier_station(m_Importer.GetBroker(), file, pierIdx, piers[pierIdx]), get_referent_direction(piers[pierIdx]));
   }
   std::sort(pier_locations.begin(), pier_locations.end(), [](const auto& a, const auto& b) {return a.first < b.first; });

   std::vector<Float64> pier_stations;
   std::vector<std::optional<Eigen::Vector2d>> referent_directions;
   for (const auto& [station, direction] : pier_locations)
   {
      pier_stations.push_back(station);
      referent_directions.push_back(direction);
   }

   for (PierIndexType pierIdx = 0; pierIdx < nPiers; pierIdx++)
   {
      bridge_desc.GetPier(pierIdx)->SetStation(pier_stations[pierIdx]);
   }

   //
   // Girders
   //

   // Locate the girders from the model geometry. The location determines the span and girder index
   auto layout = get_girder_layout(m_Importer.GetBroker(), file, pier_stations);
   locate_bearings(m_Importer.GetBroker(), file, layout);

   bool bSameNumGirdersInAllSpans = std::all_of(layout.begin(), layout.end(), [&layout](const auto& girders) {return girders.size() == layout.front().size(); });
   for (SpanIndexType spanIdx = 0; spanIdx < nSpans; spanIdx++)
   {
      if (layout[spanIdx].size() < 2)
      {
         std::_tostringstream os;
         os << _T("Span ") << LABEL_SPAN(spanIdx) << _T(" has ") << layout[spanIdx].size() << _T(" girders. At least two are required. The default number of girders is used.");
         WBFL::System::Logger::Info(os.str().c_str());
      }
   }

   if (bSameNumGirdersInAllSpans && 2 <= layout.front().size())
   {
      bridge_desc.UseSameNumberOfGirdersInAllGroups(true);
      bridge_desc.SetGirderCount(layout.front().size());
   }
   else
   {
      bridge_desc.UseSameNumberOfGirdersInAllGroups(false);
      for (SpanIndexType spanIdx = 0; spanIdx < nSpans; spanIdx++)
      {
         if (2 <= layout[spanIdx].size())
            bridge_desc.GetGirderGroup(spanIdx)->SetGirderCount(layout[spanIdx].size());
      }
   }

   // Pier orientation. The CL pier directions are needed for spacing measured along the CL piers
   auto pier_directions = set_pier_orientation(m_Importer.GetBroker(), layout, pier_stations, referent_directions, bridge_desc);

   // Girder spacing at the pier lines, measured normal to the alignment and along the CL piers
   auto [start_spacing, end_spacing] = get_beam_spacing(m_Importer.GetBroker(), layout, pier_stations);
   auto [start_spacing_along_pier, end_spacing_along_pier] = get_beam_spacing(m_Importer.GetBroker(), layout, pier_stations, &pier_directions);

   // Uniform spacing if every span has girders and all of the spacings are the same, normal to the alignment
   // or along the CL piers. Normal to the alignment first so the datum doesn't change for bridges without skew.
   bool bAllSpansHaveGirders = std::all_of(layout.begin(), layout.end(), [](const auto& girders) {return 2 <= girders.size(); });
   auto uniform_spacing = [&](const Spacing& start, const Spacing& end) -> std::optional<Float64>
      {
         std::vector<Float64> all_spacings;
         for (SpanIndexType spanIdx = 0; spanIdx < nSpans; spanIdx++)
         {
            all_spacings.insert(all_spacings.end(), start.at(spanIdx).begin(), start.at(spanIdx).end());
            all_spacings.insert(all_spacings.end(), end.at(spanIdx).begin(), end.at(spanIdx).end());
         }
         if (!bAllSpansHaveGirders || all_spacings.empty())
            return std::nullopt;
         auto [min_spacing, max_spacing] = std::minmax_element(all_spacings.begin(), all_spacings.end());
         if (WBFL::Units::ConvertToSysUnits(1.0 / 16.0, WBFL::Units::Measure::Inch) < *max_spacing - *min_spacing)
            return std::nullopt;
         return std::accumulate(all_spacings.begin(), all_spacings.end(), 0.0) / all_spacings.size();
      };

   bool bUniformSpacing = false;
   for (auto [measurement, spacing, description] : { std::make_tuple(pgsTypes::MeasurementType::NormalToItem, uniform_spacing(start_spacing, end_spacing), _T("normal to the alignment")),
                                                     std::make_tuple(pgsTypes::MeasurementType::AlongItem, uniform_spacing(start_spacing_along_pier, end_spacing_along_pier), _T("along the CL piers")) })
   {
      if (!spacing)
         continue;

      bridge_desc.SetGirderSpacingType(pgsTypes::SupportedBeamSpacing::sbsUniform);
      bridge_desc.SetGirderSpacing(*spacing);
      bridge_desc.SetMeasurementType(measurement);
      bridge_desc.SetMeasurementLocation(pgsTypes::MeasurementLocation::AtPierLine);
      bUniformSpacing = true;

      std::_tostringstream os;
      os << _T("Uniform girder spacing ") << std::fixed << std::setprecision(3) << WBFL::Units::ConvertFromSysUnits(*spacing, WBFL::Units::Measure::Feet) << _T(" ft, measured ") << description << _T(" at the pier lines");
      WBFL::System::Logger::Info(os.str().c_str());
      break;
   }

   if (!bUniformSpacing)
   {
      bridge_desc.SetGirderSpacingType(pgsTypes::SupportedBeamSpacing::sbsGeneral);
   }

   for (SpanIndexType spanIdx = 0; bUniformSpacing == false && spanIdx < nSpans; spanIdx++)
   {
      if (layout[spanIdx].size() < 2)
         continue;

      auto pSpan = bridge_desc.GetSpan(spanIdx);
      for (auto [pPier, face, spacing] : { std::make_tuple(pSpan->GetPrevPier(), pgsTypes::Ahead, &start_spacing[spanIdx]), std::make_tuple(pSpan->GetNextPier(), pgsTypes::Back, &end_spacing[spanIdx]) })
      {
         auto girder_spacing = pPier->GetGirderSpacing(face);
         girder_spacing->SetMeasurementType(pgsTypes::MeasurementType::NormalToItem);
         girder_spacing->SetMeasurementLocation(pgsTypes::MeasurementLocation::AtPierLine);
         girder_spacing->ExpandAll();
         IndexType idx = 0;
         for (auto s : *spacing)
         {
            girder_spacing->SetGirderSpacing(idx++, s);
         }
      }
   }

   SetGirderProperties(file, bridge_desc, layout);

   // Girder end distance and bearing offset normal to the CL pier
   set_end_distance_and_bearing_offset(m_Importer.GetBroker(), layout, pier_stations, pier_directions, bridge_desc);

   // Bearing data from the bearing geometry and properties
   set_bearing_data(m_Importer.GetBroker(), file, layout, bridge_desc);

   ImportSlab(file, bridge_desc, layout, pier_stations);

   pIBridgeDesc->SetBridgeDescription(bridge_desc);

   Experiment(file);

   return CIfcImporter::ImportResult::Success;
}



void CIfcBridgeImporter::SetGirderProperties(ifcopenshell::file& file, CBridgeDescription2& bridge_desc, const GirderLayout& layout)
{
   USES_CONVERSION;

   // library entry for each girder
   std::map<CGirderKey, const GirderLibraryEntry*> library_entries;
   for (const auto& girders : layout)
   {
      if (girders.size() < 2)
         continue; // girder count wasn't set for this span

      for (const auto& girder : girders)
      {
         auto beam = file.instance_by_id(girder.beam_id).as<IfcSchema::IfcBeam>();
         library_entries[girder.girder_key] = GetGirderLibraryEntry(beam);
      }
   }

   if (library_entries.empty())
      return;

   const GirderLibraryEntry* first_entry = library_entries.begin()->second;
   bool bSameGirder = std::all_of(library_entries.begin(), library_entries.end(), [first_entry](const auto& item) {return item.second == first_entry; });
   bridge_desc.UseSameGirderForEntireBridge(bSameGirder);
   if (bSameGirder)
      bridge_desc.SetGirderLibraryEntry(first_entry);

   bridge_desc.SetGirderFamilyName(first_entry->GetGirderFamilyName().c_str());

#pragma Reminder("WORKING HERE - This is assuming the first supported orientation. The IFC file doesn't have this information.")
   // should get the orientation from the girder geometry and then try to match it with the supported orientations
   // if no match, then log warning and use the default.
   auto factory = first_entry->GetBeamFactory();
   auto orientations = factory->GetSupportedGirderOrientation();
   bridge_desc.SetGirderOrientation(orientations.front());

   for (const auto& [girder_key, library_entry] : library_entries)
   {
      auto* pGirder = bridge_desc.GetGirderGroup(girder_key.groupIndex)->GetGirder(girder_key.girderIndex);
      if (!bSameGirder)
         pGirder->SetGirderLibraryEntry(library_entry);

      const auto& girder = layout[girder_key.groupIndex][girder_key.girderIndex];
      auto beam = file.instance_by_id(girder.beam_id).as<IfcSchema::IfcBeam>();

      auto fci = GetMeasureProperty<IfcSchema::IfcPressureMeasure>(CIfcImporter::GetUnits(), beam, "Pset_PrecastConcreteElementGeneral", "ReleaseStrength");
      if (fci)
      {
         pGirder->GetSegment(0)->Material.Concrete.Fci = *fci;
      }
      else
      {
         WBFL::System::Logger::Info(_T("ReleaseStrength property in Pset_PrecastConcreteElementGeneral property set not found"));
      }

      auto material = GetMaterial<IfcSchema>(beam);
      if (material)
      {
         auto fc = GetMaterialMeasureProperty<IfcSchema::IfcPressureMeasure>(CIfcImporter::GetUnits(), material, "Pset_MaterialConcrete", "CompressiveStrength");
         if (fc)
         {
            pGirder->GetSegment(0)->Material.Concrete.Fc = *fc;
         }
         else
         {
            WBFL::System::Logger::Info(_T("CompressiveStrength property in Pset_MaterialConcrete property set not found"));
         }
      }
      else
      {
         WBFL::System::Logger::Info(_T("Materials are not associated with the beam"));
      }
   }
}

bool CIfcBridgeImporter::IsValidBridge(ifcopenshell::file& file, IfcSchema::IfcBridge bridge)
{
   // must be a girder bridge
   if (bridge.PredefinedType().value_or(IfcSchema::IfcBridgeTypeEnum::IfcBridgeType_NOTDEFINED) != IfcSchema::IfcBridgeTypeEnum::IfcBridgeType_GIRDER)
      return false;

   // This check could be far less strict if we can assume the user provided us a PSG bridge.
   // If we can go from the girder line geometry and the girder name, mapped to a library entry,
   // that might be enough to actually do some work
   //if (!HasValidGirders(file, bridge))
   //   return false;

   return true;
}

bool CIfcBridgeImporter::HasValidGirders(ifcopenshell::file& file, IfcSchema::IfcBridge bridge)
{
   if (HasValidGirdersByTPF(file, bridge))
      return true;

   if (HasValidGirdersByOther(file, bridge))
      return true;

   WBFL::System::Logger::Info(_T("One or more beams in the superstructure could not be identified as precast, prestressed concrete."));

   return false;
}

bool CIfcBridgeImporter::HasValidGirdersByTPF(ifcopenshell::file& file, IfcSchema::IfcBridge bridge)
{
   auto superstructure = GetBridgePart(file, IfcSchema::IfcBridgePartTypeEnum::IfcBridgePartType_SUPERSTRUCTURE);

   auto beams = file.instances_by_type<IfcSchema::IfcBeam>();
   bool valid_beams = true;
   for (auto& beam : beams)
   {
      // beam must be contained in the spatial structure of the superstructure
      auto related_elements = beam.ContainedInStructure();
      if (!related_elements.empty())
      {
         for (auto& related_element : related_elements)
         {
            if (related_element.RelatingStructure() != superstructure)
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
            if (!HasClassification<IfcSchema>(beam, "GirderPrestressedConcrete"))
            {
               valid_beams = false;
               break;
            }
         }
      }
   }

   return valid_beams;
}

bool CIfcBridgeImporter::HasValidGirdersByOther(ifcopenshell::file& file, IfcSchema::IfcBridge bridge)
{
   // This function attempts to determine if all the beams in the superstructure are precast concrete.
   // The basic idea is that the beams are IfcElementAssembly and they are factory assembled girders.
   // A factory assembled girder alone is not enough to claim the girders are precast concrete.
   //
   // This function may be a bad idea... keep it for now, but be skeptical
   auto superstructure = GetBridgePart(file, IfcSchema::IfcBridgePartTypeEnum::IfcBridgePartType_SUPERSTRUCTURE);

   auto element_assemblies = file.instances_by_type<IfcSchema::IfcElementAssembly>();
   bool valid_beams = true;
   for (auto& element_assembly : element_assemblies)
   {
      // beam must be contained in the spatial structure of the superstructure
      auto related_elements = element_assembly.ContainedInStructure();
      if (!related_elements.empty())
      {
         for (auto& related_element : related_elements)
         {
            if (related_element.RelatingStructure() != superstructure)
               continue;
         }
      }
      else
      {
         continue; // not in a spatial structure
      }

      auto assembly_place = element_assembly.AssemblyPlace();
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

IfcSchema::IfcBridge CIfcBridgeImporter::GetBridge(ifcopenshell::file& file)
{
   USES_CONVERSION;

   auto bridges = file.instances_by_type<IfcSchema::IfcBridge>();

   if (1 <= bridges.size())
   {
      std::vector<IfcSchema::IfcBridge> valid_bridges;
      for (auto& bridge : bridges)
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
               auto strLabel = (bridge.Name() ? *(bridge.Name()) : bridge.Description() ? *(bridge.Description()) : "Unnamed");
               os << strLabel << std::endl;
            }

            if (m_Importer.IsInteractive())
            {
               result = AfxChoose(_T("Select Bridge"), _T("Select bridge to import"), A2T(os.str().c_str()), 0, TRUE);
               if (result < 0)
                  return {};
            }
            else
            {
               WBFL::System::Logger::Info(_T("IFC model contains more than one bridge. Importing the first bridge."));
            }
         }

         return valid_bridges[result];
      }
   }

   return {};
}

void CIfcBridgeImporter::ImportSlab(ifcopenshell::file& file, CBridgeDescription2& bridge_desc, const GirderLayout& layout, const std::vector<Float64>& pier_stations)
{
   auto slabs = file.instances_by_type<IfcSchema::IfcSlab>();
   auto it = std::find_if(slabs.begin(), slabs.end(), [](const auto& slab) {return slab.PredefinedType() == IfcSchema::IfcSlabTypeEnum::IfcSlabType_FLOOR; });
   if (it == slabs.end())
      return; // no slabs

   auto deck = get_deck_mesh(m_Importer.GetBroker(), file);
   if (!deck)
      return; // keep the default deck

   auto* pDeck = bridge_desc.GetDeckDescription();

   // Stay-in-place deck panels are not supported yet
   pDeck->SetDeckType(pgsTypes::sdtCompositeCIP);

   // plan view shape
   auto deck_points = get_deck_edge_points(m_Importer.GetBroker(), *deck);
   if (!deck_points.empty())
      pDeck->DeckEdgePoints = deck_points;

   // cross section, measured with the haunches, which may be separate from the deck slab
   auto deck_concrete = get_deck_concrete_mesh(m_Importer.GetBroker(), file);
   auto section = analyze_deck_section(m_Importer.GetBroker(), deck_concrete ? *deck_concrete : *deck, layout, pier_stations, bridge_desc);

   auto log_length = [](LPCTSTR name, Float64 value)
      {
         std::_tostringstream os;
         os << name << _T(" ") << std::fixed << std::setprecision(3) << WBFL::Units::ConvertFromSysUnits(value, WBFL::Units::Measure::Inch) << _T(" in (from deck geometry)");
         WBFL::System::Logger::Info(os.str().c_str());
      };
   auto log_default = [](LPCTSTR name)
      {
         std::_tostringstream os;
         os << name << _T(" could not be measured from the deck geometry. The default is used.");
         WBFL::System::Logger::Info(os.str().c_str());
      };

   if (section.gross_depth)
   {
      pDeck->GrossDepth = *section.gross_depth;
      log_length(_T("Deck gross depth"), *section.gross_depth);
   }
   else
   {
      log_default(_T("Deck gross depth"));
   }

   for (auto side : { pgsTypes::stLeft, pgsTypes::stRight })
   {
      LPCTSTR strSide = (side == pgsTypes::stLeft ? _T("Left") : _T("Right"));
      std::_tstring edge_name = std::_tstring(strSide) + _T(" overhang edge depth");
      if (section.edge_depth[side])
      {
         pDeck->OverhangEdgeDepth[side] = *section.edge_depth[side];
         log_length(edge_name.c_str(), *section.edge_depth[side]);
      }
      else
      {
         log_default(edge_name.c_str());
      }

      std::_tstring taper_name = std::_tstring(strSide) + _T(" overhang taper");
      if (section.overhang_taper[side])
      {
         pDeck->OverhangTaper[side] = *section.overhang_taper[side];
         static const std::map<pgsTypes::DeckOverhangTaper, LPCTSTR> taper_names{
            {pgsTypes::dotNone, _T("none")},
            {pgsTypes::dotTopTopFlange, _T("taper to top of girder top flange")},
            {pgsTypes::dotBottomTopFlange, _T("taper to bottom of girder top flange")} };
         std::_tostringstream os;
         os << taper_name << _T(": ") << taper_names.at(*section.overhang_taper[side]) << _T(" (from deck geometry)");
         WBFL::System::Logger::Info(os.str().c_str());
      }
      else
      {
         log_default(taper_name.c_str());
      }
   }

   if (section.haunch_shape)
   {
      pDeck->HaunchShape = *section.haunch_shape;
      WBFL::System::Logger::Info(*section.haunch_shape == pgsTypes::hsSquare ? _T("Haunch shape: square (from deck geometry)") : _T("Haunch shape: filleted (from deck geometry)"));
      if (section.fillet)
      {
         bridge_desc.SetFillet(*section.fillet);
         log_length(_T("Fillet"), *section.fillet);
      }
   }
   else
   {
      WBFL::System::Logger::Info(_T("The haunches are not part of the deck geometry. The default haunch shape and fillet are used."));
   }

   // Slab offset ("A" dimension) at CL bearings
   std::vector<Float64> all_slab_offsets;
   for (const auto& girders : section.slab_offset)
   {
      for (const auto& ends : girders)
      {
         for (const auto& value : ends)
         {
            if (value)
               all_slab_offsets.push_back(*value);
         }
      }
   }

   if (all_slab_offsets.empty())
   {
      log_default(_T("Slab offset"));
      return;
   }

   bridge_desc.SetHaunchInputDepthType(pgsTypes::hidACamber);

   if (section.slab_offset_extrapolated)
      WBFL::System::Logger::Info(_T("The deck doesn't extend over some bearings (e.g. integral abutments). The slab offset at those bearings is extrapolated from where the deck is over the girder."));

   const Float64 tolerance = WBFL::Units::ConvertToSysUnits(1.0 / 16.0, WBFL::Units::Measure::Inch);
   auto spread = [](const std::vector<Float64>& values) {auto [lo, hi] = std::minmax_element(values.begin(), values.end()); return *hi - *lo; };
   auto average = [](const std::vector<Float64>& values) {return std::accumulate(values.begin(), values.end(), 0.0) / values.size(); };

   // slab offsets on each bearing line (pier face)
   std::map<std::pair<PierIndexType, pgsTypes::PierFaceType>, std::vector<Float64>> bearing_lines;
   for (SpanIndexType spanIdx = 0; spanIdx < section.slab_offset.size(); spanIdx++)
   {
      for (const auto& ends : section.slab_offset[spanIdx])
      {
         if (ends[pgsTypes::metStart]) bearing_lines[{spanIdx, pgsTypes::Ahead}].push_back(*ends[pgsTypes::metStart]);
         if (ends[pgsTypes::metEnd]) bearing_lines[{spanIdx + 1, pgsTypes::Back}].push_back(*ends[pgsTypes::metEnd]);
      }
   }

   if (spread(all_slab_offsets) <= tolerance)
   {
      bridge_desc.SetSlabOffsetType(pgsTypes::sotBridge);
      bridge_desc.SetSlabOffset(average(all_slab_offsets));
      log_length(_T("Slab offset for the bridge"), average(all_slab_offsets));
   }
   else if (std::all_of(bearing_lines.begin(), bearing_lines.end(), [&](const auto& item) {return spread(item.second) <= tolerance; }))
   {
      bridge_desc.SetSlabOffsetType(pgsTypes::sotBearingLine);
      for (const auto& [key, values] : bearing_lines)
      {
         auto [pierIdx, face] = key;
         bridge_desc.GetPier(pierIdx)->SetSlabOffset(face, average(values));
         std::_tostringstream os;
         os << _T("Slab offset at ") << LABEL_PIER(pierIdx) << (face == pgsTypes::Back ? _T(" back") : _T(" ahead"));
         log_length(os.str().c_str(), average(values));
      }
   }
   else
   {
      bridge_desc.SetSlabOffsetType(pgsTypes::sotSegment);
      for (SpanIndexType spanIdx = 0; spanIdx < section.slab_offset.size(); spanIdx++)
      {
         for (GirderIndexType gdrIdx = 0; gdrIdx < section.slab_offset[spanIdx].size(); gdrIdx++)
         {
            const auto& ends = section.slab_offset[spanIdx][gdrIdx];
            if (layout[spanIdx].size() < 2)
               continue; // girders weren't imported for this span

            auto* pSegment = bridge_desc.GetGirderGroup(spanIdx)->GetGirder(gdrIdx)->GetSegment(0);
            for (auto end : { pgsTypes::metStart, pgsTypes::metEnd })
            {
               if (ends[end])
                  pSegment->SetSlabOffset(end, *ends[end]);
            }
         }
      }
      WBFL::System::Logger::Info(_T("Slab offset varies by girder. The slab offset is set for each segment (from deck and girder geometry)."));
   }
}

namespace
{
   // upper case letters and digits only, so "WF66G", "wf 66 g" and "WF-66G" compare equal
   std::string normalize_girder_name(const std::string& name)
   {
      std::string result;
      for (auto c : name)
      {
         if (std::isalnum(static_cast<unsigned char>(c)))
            result.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
      }
      return result;
   }

   // similarity of two normalized names in the range [0,1], based on the longest common substring
   double girder_name_similarity(const std::string& a, const std::string& b)
   {
      if (a.empty() || b.empty())
         return 0.0;

      std::vector<size_t> prev(b.size() + 1, 0), curr(b.size() + 1, 0);
      size_t longest = 0;
      for (size_t i = 1; i <= a.size(); i++)
      {
         for (size_t j = 1; j <= b.size(); j++)
         {
            curr[j] = (a[i - 1] == b[j - 1]) ? prev[j - 1] + 1 : 0;
            longest = std::max(longest, curr[j]);
         }
         std::swap(prev, curr);
      }
      return (double)longest / (double)std::max(a.size(), b.size());
   }

   // true if a property or classification name suggests it holds the girder type (e.g. "Type", "2_Type", "ShapeName")
   bool is_girder_type_property(const std::string& name)
   {
      auto n = normalize_girder_name(name);
      return n.ends_with("TYPE") || n.find("SHAPE") != std::string::npos;
   }
}

std::vector<std::string> CIfcBridgeImporter::GetGirderTypeNames(IfcSchema::IfcBeam beam)
{
   // Collect the names that might identify the girder type, most authoritative first.
   // Models that don't follow the TPF/usBridge conventions put the girder type in all sorts of places
   std::vector<std::string> names;
   auto add = [&names](const std::string& name)
      {
         // skip empty names and generic words that say nothing about the girder type
         static const std::set<std::string> generic{ "", "NULL", "BEAM", "GIRDER", "NOTDEFINED", "USERDEFINED" };
         if (generic.find(normalize_girder_name(name)) == generic.end() && std::find(names.begin(), names.end(), name) == names.end())
            names.push_back(name);
      };

   auto type = GetType<IfcSchema::IfcBeamType>(beam);
   if (type && type.Name())
      add(*type.Name());

   if (beam.ObjectType())
      add(*beam.ObjectType());

   // text valued properties with names like Type, 2_Type, or ShapeName
   for (auto& rel : beam.IsDefinedBy())
   {
      auto pset = rel.RelatingPropertyDefinition().as<IfcSchema::IfcPropertySet>();
      if (!pset)
         continue;

      for (auto& property : pset.HasProperties())
      {
         auto value = property.as<IfcSchema::IfcPropertySingleValue>();
         if (!value || !is_girder_type_property(value.Name()) || !value.NominalValue())
            continue;

         if (auto label = value.NominalValue().as<IfcSchema::IfcLabel>())
            add(label);
         else if (auto text = value.NominalValue().as<IfcSchema::IfcText>())
            add(text);
      }
   }

   // classification references (e.g. "Beam, PPC, BTB45")
   for (auto& rel : beam.HasAssociations())
   {
      auto rel_classification = rel.as<IfcSchema::IfcRelAssociatesClassification>();
      auto reference = rel_classification ? rel_classification.RelatingClassification().as<IfcSchema::IfcClassificationReference>() : IfcSchema::IfcClassificationReference{};
      // usBridge classifications (e.g. usBridge_GirderPrecastConcrete) identify the kind of element, not the girder type
      if (reference && reference.Name() && !reference.Identification().value_or("").starts_with("usBridge_"))
         add(*reference.Name());
   }

   return names;
}

const GirderLibraryEntry* CIfcBridgeImporter::GetGirderLibraryEntry(IfcSchema::IfcBeam beam)
{
   USES_CONVERSION;

   auto girder_type_names = GetGirderTypeNames(beam);
   auto found = m_GirderMatches.find(girder_type_names);
   if (found != m_GirderMatches.end())
      return found->second;

   GET_IFACE2(m_Importer.GetBroker(), ILibrary, pLibrary);
   GET_IFACE2(m_Importer.GetBroker(), ILibraryNames, pLibNames);
   std::vector<std::_tstring> library_names;
   pLibNames->EnumGirderNames(&library_names); // all girder families

   // Exact match on any of the names (ignoring case, spaces, and punctuation).
   // Otherwise, use the library girder with the most similar name.
   // Ultimately, girders that aren't in the library should be created, or the user
   // prompted to pick or create one. See devdocs/IfcImportPlan.md
   const GirderLibraryEntry* girder_library_entry = nullptr;
   std::_tstring best_name;
   double best_similarity = -1.0;
   for (const auto& girder_type_name : girder_type_names)
   {
      auto normalized_type_name = normalize_girder_name(girder_type_name);
      for (const auto& library_name : library_names)
      {
         double similarity = girder_name_similarity(normalized_type_name, normalize_girder_name(T2A(library_name.c_str())));
         if (best_similarity < similarity)
         {
            best_similarity = similarity;
            best_name = library_name;
         }
      }
   }

   // Names this dissimilar aren't really a match. Fall back to the first I-Beam, or the first girder of any kind
   const double min_similarity = 0.5;
   bool bFallback = (best_similarity < min_similarity);
   if (bFallback)
   {
      std::vector<std::_tstring> ibeam_names;
      pLibNames->EnumGirderNames(_T("I-Beam"), &ibeam_names);
      best_name = !ibeam_names.empty() ? ibeam_names.front() : (!library_names.empty() ? library_names.front() : std::_tstring());
   }

   if (!best_name.empty())
      girder_library_entry = pLibrary->GetGirderEntry(best_name.c_str());

   if (best_similarity < 1.0)
   {
      std::ostringstream os;
      os << "Girder type of " << beam.Name().value_or("unnamed beam") << " ";
      if (girder_type_names.empty())
      {
         os << "could not be determined";
      }
      else
      {
         os << "(";
         for (const auto& name : girder_type_names)
            os << "\"" << name << "\"" << (&name == &girder_type_names.back() ? "" : ", ");
         os << ") was not found in the library";
      }
      if (bFallback)
         os << ". No similar girder in the library. Using " << (best_name.empty() ? "(none)" : T2A(best_name.c_str()));
      else
         os << ". Using best match " << T2A(best_name.c_str()) << " (similarity " << std::fixed << std::setprecision(2) << best_similarity << ")";

      WBFL::System::Logger::Info(os.str().c_str());
   }

   m_GirderMatches.emplace(girder_type_names, girder_library_entry);
   return girder_library_entry;
}

bool CIfcBridgeImporter::DeriveAlignmentFromDeck(ifcopenshell::file& file)
{
   return create_alignment_from_deck(m_Importer.GetBroker(), file);
}

void CIfcBridgeImporter::Experiment(ifcopenshell::file& file)
{
   //get_beam_spacing(file);
}

