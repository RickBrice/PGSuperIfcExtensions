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
#include "IfcTargetReader.h"

#include "BeamSpacing.h"
#include "DeckSlab.h"
#include "TrafficBarriers.h"
#include "RoadwaySection.h"
#include "Piers.h"
#include "Strands.h"

#include <chrono>
#include <numeric>

#include <MFCTools\Prompts.h>
#include <EAF/EAFProgress.h>
#include <boost/range/combine.hpp>
#include <psgLib/BridgeDescription2.h>
#include <psgLib/GirderLabel.h>
#include <psgLib/GirderLibraryEntry.h>
#include <psgLib/LibraryManager.h>
#include <psgLib/BeamFamilyManager.h>
#include <IFace/DocumentType.h>
#include <PGSuperCatCom.h>
#include <PGSpliceCatCom.h>



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

   // The number of spans is from the modeled substructure. A number of spans property is only a cross check;
   // when it disagrees, the model is used and the conflict is reported (devdocs/MappingTablesDesign.md, G4)
   SpanIndexType nSpans = nPiers - 1;
   const auto& reader = CIfcImporter::GetTargetReader();
   if (auto reading = reader.Read("bridge.number_of_spans", bridge))
   {
      auto nSpansProperty = std::get<Int64>(reading->value);
      if (nSpansProperty != (Int64)nSpans)
      {
         std::ostringstream os;
         os << "The number of spans is " << nSpansProperty << " (" << reading->Source() << "), but the model has " << nPiers
            << " abutments and piers, so " << nSpans << " spans. The modeled substructure is used.";
         WBFL::System::Logger::Error(os.str());
      }
   }
   else
   {
      reader.ReportNotFound("bridge.number_of_spans", bridge, "the bridge");
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

   // Prestressing strands, after the library entries are set
   import_strands(m_Importer.GetBroker(), file, layout, bridge_desc);

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
         library_entries[girder.girder_key] = GetGirderLibraryEntry(beam, girder);
      }
   }

   if (library_entries.empty())
      return;

   const GirderLibraryEntry* first_entry = library_entries.begin()->second;
   bool bSameGirder = std::all_of(library_entries.begin(), library_entries.end(), [first_entry](const auto& item) {return item.second == first_entry; });
   bridge_desc.UseSameGirderForEntireBridge(bSameGirder);
   if (bSameGirder)
   {
      // the entry and its name: PGSuper finds the entry by name when the project is loaded (SetGirderLibraryEntry
      // doesn't set the name, unlike CSplicedGirderData's). The entry first, since setting the name copies the entry's
      // shear and longitudinal reinforcement to the segments
      bridge_desc.SetGirderLibraryEntry(first_entry);
      bridge_desc.SetGirderName(first_entry->GetName().c_str());
   }

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
      auto beam_name = beam.Name().value_or("unnamed beam");

      const auto& reader = CIfcImporter::GetTargetReader();
      for (auto [target, label, pStrength] : { std::make_tuple("girder.fci", "f'ci", &pGirder->GetSegment(0)->Material.Concrete.Fci),
                                               std::make_tuple("girder.fc", "f'c", &pGirder->GetSegment(0)->Material.Concrete.Fc) })
      {
         if (auto reading = reader.Read(target, beam))
         {
            *pStrength = std::get<Float64>(reading->value);

            std::ostringstream os;
            os << label << " of " << beam_name << " = " << std::fixed << std::setprecision(3) << WBFL::Units::ConvertFromSysUnits(*pStrength, WBFL::Units::Measure::KSI) << " ksi, " << reading->Source();
            WBFL::System::Logger::Info(os.str());
         }
         else
         {
            reader.ReportNotFound(target, beam, beam_name);
         }
      }
   }

   // The girders of a girder type group share a library entry (PGSuper checks it). With different entries, each run of
   // adjacent girders with the same entry is a group
   if (!bSameGirder)
   {
      for (GroupIndexType grpIdx = 0; grpIdx < bridge_desc.GetGirderGroupCount(); grpIdx++)
      {
         auto* pGroup = bridge_desc.GetGirderGroup(grpIdx);
         pGroup->ExpandAll();
         GirderIndexType nGirders = pGroup->GetGirderCount();
         for (GirderIndexType first = 0; first < nGirders; )
         {
            std::_tstring name(pGroup->GetGirder(first)->GetGirderName());
            GirderIndexType last = first;
            while (last + 1 < nGirders && name == pGroup->GetGirder(last + 1)->GetGirderName())
               last++;
            if (first < last)
               pGroup->Join(first, last, first);
            first = last + 1;
         }
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
         const auto& reader = CIfcImporter::GetTargetReader();
         auto assembly_place = reader.Read("girder.assembly_place", beam);
         auto casting_method = reader.Read("girder.casting_method", beam);

         // level 1 check - if there is an assembly_place and casting_method, use that for the determination
         // otherwise use the more specific usBridge classification
         if (assembly_place && casting_method)
         {
            bool is_factory_assembled = (std::get<std::string>(assembly_place->value) == std::string("FACTORY"));
            bool is_precast = (std::get<std::string>(casting_method->value) == std::string("PRECAST"));

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
   if (get_slab_id(file) == 0)
      return; // no deck

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

   // the profile (when the alignment has none) and the roadway section, from the top of the deck between the barriers,
   // then the traffic barriers, which PGSuper places at the edges of the deck turned to the roadway's cross slope
   auto barriers = get_barrier_meshes(m_Importer.GetBroker(), file);
   import_profile(m_Importer.GetBroker(), *deck, barriers, pier_stations, !m_Importer.HasProfile());
   import_roadway_section(m_Importer.GetBroker(), *deck, barriers, pier_stations);
   import_traffic_barriers(m_Importer.GetBroker(), *deck, barriers, pier_stations, bridge_desc);

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

   // A deck depth property is a cross check: geometry wins and a conflict is reported (D1).
   // The property is used when the depth can't be measured from the geometry
   auto deck_slab = file.instance_by_id(get_slab_id(file)).as<IfcSchema::IfcObject>();
   const auto& reader = CIfcImporter::GetTargetReader();
   auto depth_reading = reader.Read("deck.gross_depth", deck_slab);
   if (section.gross_depth)
   {
      pDeck->GrossDepth = *section.gross_depth;
      log_length(_T("Deck gross depth"), *section.gross_depth);

      if (depth_reading)
      {
         Float64 property_depth = std::get<Float64>(depth_reading->value);
         const Float64 tolerance = WBFL::Units::ConvertToSysUnits(1.0 / 8.0, WBFL::Units::Measure::Inch);
         if (!IsEqual(property_depth, *section.gross_depth, tolerance))
         {
            std::ostringstream os;
            os << std::fixed << std::setprecision(3) << "The deck gross depth is " << WBFL::Units::ConvertFromSysUnits(property_depth, WBFL::Units::Measure::Inch)
               << " in (" << depth_reading->Source() << "), but measures " << WBFL::Units::ConvertFromSysUnits(*section.gross_depth, WBFL::Units::Measure::Inch)
               << " in from the deck geometry. The geometry is used.";
            WBFL::System::Logger::Warning(os.str());
         }
      }
   }
   else if (depth_reading)
   {
      pDeck->GrossDepth = std::get<Float64>(depth_reading->value);
      std::ostringstream os;
      os << std::fixed << std::setprecision(3) << "Deck gross depth " << WBFL::Units::ConvertFromSysUnits(pDeck->GrossDepth, WBFL::Units::Measure::Inch)
         << " in (" << depth_reading->Source() << "). It couldn't be measured from the deck geometry.";
      WBFL::System::Logger::Info(os.str());
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

}

std::vector<std::string> CIfcBridgeImporter::GetGirderTypeNames(IfcSchema::IfcBeam beam)
{
   // The names that identify the girder type, from the locations of girder.type_names in the mapping table, in table order
   std::vector<std::string> names;
   for (const auto& reading : CIfcImporter::GetTargetReader().ReadAll("girder.type_names", beam))
   {
      // skip generic words that say nothing about the girder type
      static const std::set<std::string> generic{ "", "NULL", "BEAM", "GIRDER", "NOTDEFINED", "USERDEFINED" };
      const auto& name = std::get<std::string>(reading.value);
      if (generic.find(normalize_girder_name(name)) == generic.end())
         names.push_back(name);
   }
   return names;
}

const GirderLibraryEntry* CIfcBridgeImporter::MatchLibraryEntry(const std::string& beam_name, const std::vector<std::string>& girder_type_names, const SectionMeasures& section)
{
   USES_CONVERSION;
   GET_IFACE2_NOCHECK(m_Importer.GetBroker(), ILibrary, pLibrary); // not used when the library sections are cached and nothing matches
   GET_IFACE2(m_Importer.GetBroker(), ILibraryNames, pLibNames);
   std::vector<std::_tstring> library_names;
   pLibNames->EnumGirderNames(&library_names); // all girder families

   std::set<std::string> normalized_type_names;
   for (const auto& name : girder_type_names)
      normalized_type_names.insert(normalize_girder_name(name));

   // the library entries whose sections match, and the closest entry
   std::vector<std::pair<std::_tstring, double>> matches; // name, difference
   std::_tstring closest_name;
   double closest_difference = DBL_MAX;
   std::vector<std::string> closest_differences;
   std::_tstring named_entry; // a library entry with one of the girder type names
   std::vector<std::string> named_differences;
   std::vector<std::pair<double, std::_tstring>> ranked; // difference, name
   std::vector<std::string> template_entries;
   for (const auto& library_name : library_names)
   {
      if (!IsUsableLibraryEntry(library_name))
      {
         template_entries.push_back(T2A(library_name.c_str()));
         continue;
      }

      auto found = m_LibrarySections.find(library_name);
      if (found == m_LibrarySections.end())
      {
         std::optional<SectionMeasures> measures;
         if (auto outline = GetLibraryEntryOutline(pLibrary->GetGirderEntry(library_name.c_str())))
            measures = MeasureSection(*outline);
         found = m_LibrarySections.emplace(library_name, measures).first;
      }
      if (!found->second)
         continue;

      auto differences = m_pMatchComparison->Compare(section, *found->second);
      double difference = m_pMatchComparison->Difference(section, *found->second);
      ranked.emplace_back(difference, library_name);
      if (differences.empty())
         matches.emplace_back(library_name, difference);
      if (difference < closest_difference)
      {
         closest_difference = difference;
         closest_name = library_name;
         closest_differences = differences;
      }
      if (normalized_type_names.count(normalize_girder_name(T2A(library_name.c_str()))))
      {
         named_entry = library_name;
         named_differences = differences;
      }
   }

   if (!m_bLibrarySectionsReported)
   {
      auto n = std::count_if(m_LibrarySections.begin(), m_LibrarySections.end(), [](const auto& item) {return item.second.has_value(); });
      std::ostringstream os;
      os << n << " of " << library_names.size() - template_entries.size() << " library girders have a section to compare with.";
      if (!template_entries.empty())
      {
         os << " Not compared, since they came with the import template:";
         for (const auto& name : template_entries)
            os << " " << name << (&name == &template_entries.back() ? "." : ",");
      }
      std::vector<std::string> without;
      for (const auto& [name, measures] : m_LibrarySections)
         if (!measures)
            without.push_back(T2A(name.c_str()));
      if (!without.empty())
      {
         os << " Without: ";
         for (const auto& name : without)
            os << name << (&name == &without.back() ? "" : ", ");
      }
      WBFL::System::Logger::Info(os.str());
      m_bLibrarySectionsReported = true;
   }

   std::ostringstream os;
   os << "Section of " << beam_name << ": " << DescribeSection(section) << ". ";
   const GirderLibraryEntry* entry = nullptr;
   if (!named_entry.empty() && named_differences.empty())
   {
      // the name is confirmed by the section
      os << "Library girder " << T2A(named_entry.c_str()) << " has the girder type name and the same section.";
      entry = pLibrary->GetGirderEntry(named_entry.c_str());
   }
   else if (!matches.empty())
   {
      std::sort(matches.begin(), matches.end(), [](const auto& a, const auto& b) {return a.second < b.second; });
      os << "Library girder " << T2A(matches.front().first.c_str()) << " has the same section";
      if (1 < matches.size())
         os << " (" << matches.size() << " library girders do, this one is the closest)";
      if (!named_entry.empty())
      {
         os << ". The library girder with the girder type name, " << T2A(named_entry.c_str()) << ", has a different section (";
         for (const auto& d : named_differences)
            os << d << (&d == &named_differences.back() ? "" : ", ");
         os << ")";
      }
      os << ".";
      entry = pLibrary->GetGirderEntry(matches.front().first.c_str());
   }
   else
   {
      os << "No library girder has this section.";
      if (!closest_name.empty())
      {
         os << " The closest is " << T2A(closest_name.c_str()) << " (";
         for (const auto& d : closest_differences)
            os << d << (&d == &closest_differences.back() ? "" : ", ");
         os << ").";
      }

      // G1 step 4: fit the dimensions of every beam factory (a family has several, e.g. WF, NU, and AASHTO girders are
      // I-beams), starting from the factory's closest library entry and from its default dimensions (a fit can settle
      // near its start: a 58 in entry fitted to a 31 in girder may not find the fit the defaults do). The best fit makes a
      // new library entry with a generated strand grid; entries aren't cloned (the user, 2026-10-20). Every installed
      // beam factory is tried, so no configuration is assumed
      std::sort(ranked.begin(), ranked.end());
      struct Candidate
      {
         std::_tstring entry; // the library entry to clone, or empty for the factory's defaults
         std::shared_ptr<PGS::Beams::BeamFactory> factory;
         GirderDimensions start;
      };
      std::vector<Candidate> candidates;
      std::set<std::_tstring> factory_names;
      std::map<std::_tstring, std::_tstring> factory_entries; // factory name -> its closest library entry
      for (const auto& [difference, name] : ranked)
      {
         const auto* pEntry = pLibrary->GetGirderEntry(name.c_str());
         auto factory = pEntry->GetBeamFactory();
         if (factory && factory_names.insert(factory->GetName()).second)
         {
            candidates.push_back({ name, factory, pEntry->GetDimensions() }); // the entry's dimensions are a start for the fit
            factory_entries[factory->GetName()] = name;
         }
      }
      GET_IFACE2(m_Importer.GetBroker(), IDocumentType, pDocType);
      for (const auto& family_name : PGS::Library::BeamFamilyManager::GetBeamFamilyNames(pDocType->IsPGSpliceDocument() ? CATID_PGSpliceBeamFamily : CATID_PGSuperBeamFamily))
      {
         auto family = PGS::Library::BeamFamilyManager::GetBeamFamily(family_name);
         if (!family)
            continue;
         for (const auto& factory_name : family->GetFactoryNames())
         {
            auto factory = family->CreateFactory(factory_name);
            if (!factory)
               continue;
            auto entry = factory_entries.find(factory->GetName());
            if (entry == factory_entries.end() && !factory_names.insert(factory->GetName()).second)
               continue; // already a candidate from another family
            GirderDimensions start;
            const auto& names = factory->GetDimensionNames();
            const auto& defaults = factory->GetDefaultDimensions();
            for (size_t i = 0; i < names.size() && i < defaults.size(); i++)
               start.emplace_back(names[i], defaults[i]);
            if (entry != factory_entries.end())
               factory_entries.erase(entry); // its defaults once
            candidates.push_back({ std::_tstring(), factory, start });
         }
      }

      std::optional<DimensionFit> best_fit;
      const Candidate* best_candidate = nullptr;
      std::string how; // how the fit was found, for the log
      auto fit_start_time = std::chrono::steady_clock::now();
      // An accepted fit is better than one that isn't, whatever it started from (the new entry is made by the factory
      // either way). Otherwise the smaller difference is better
      auto rank = [](const Candidate&, const DimensionFit& fit)
         {
            return fit.differences.empty() ? 0 : 2;
         };
      std::vector<std::pair<const Candidate*, DimensionFit>> fits; // every fit considered, for the ties
      auto consider = [&](const Candidate& candidate, const std::optional<DimensionFit>& fit)
         {
            if (fit)
               fits.emplace_back(&candidate, *fit);
            if (fit && (!best_fit || rank(candidate, *fit) < rank(*best_candidate, *best_fit) ||
               (rank(candidate, *fit) == rank(*best_candidate, *best_fit) && fit->difference < best_fit->difference)))
            {
               best_fit = fit;
               best_candidate = &candidate;
            }
         };

      // First the library entries whose top and bottom flanges match the model's: girders of a series (e.g. the WSDOT
      // WF girders) have the same flanges, so only the height differs
      std::vector<Candidate> flange_candidates;
      for (const auto& [difference, name] : ranked)
      {
         const auto& measures = m_LibrarySections[name];
         if (measures && FlangesMatch(section, *measures, SectionTolerances()))
         {
            const auto* pEntry = pLibrary->GetGirderEntry(name.c_str());
            flange_candidates.push_back({ name, pEntry->GetBeamFactory(), pEntry->GetDimensions() });
         }
      }
      // The model's girder has a constant depth (the section is measured over its prismatic middle), so a fit must
      // make a girder of constant depth (MakeConstantDepth: e.g. a tapered I-beam fitted at one end has the other free)
      std::set<std::string> varying_depth; // factories that can't make the section as a girder of constant depth
      auto prismatic = [&varying_depth](const std::shared_ptr<PGS::Beams::BeamFactory>& factory, const std::optional<DimensionFit>& fit) -> std::optional<DimensionFit>
         {
            if (!fit)
               return fit;
            auto result = MakeConstantDepth(factory, *fit);
            if (!result)
               varying_depth.insert(std::string(CT2A(factory->GetName().c_str())));
            return result;
         };

      for (const auto& candidate : flange_candidates)
      {
         auto fit = prismatic(candidate.factory, m_pHeightFitter->Fit(candidate.factory, candidate.start, section, *m_pFitComparison));
         if (fit && fit->differences.empty())
            consider(candidate, fit);
      }
      if (best_fit)
      {
         how = "the flanges match, only the height is fitted";
      }
      else
      {
         // Otherwise all of the dimensions of every beam factory: a short fit from each start, then a full fit of the
         // five most promising factories (each from its better start). Fitting takes a while, so the progress window
         // says which fit is running
         GET_IFACE2(m_Importer.GetBroker(), IEAFProgress, pProgress);
         std::vector<std::pair<double, const Candidate*>> screened;
         size_t count = 0;
         for (const auto& candidate : candidates)
         {
            std::_tostringstream progress;
            progress << _T("Fitting the section of ") << A2T(beam_name.c_str()) << _T(": screening ") << ++count << _T(" of ") << candidates.size() << _T(" (") << candidate.factory->GetName() << _T(")");
            pProgress->UpdateMessage(progress.str().c_str());
            if (auto fit = prismatic(candidate.factory, m_pScreenFitter->Fit(candidate.factory, candidate.start, section, *m_pFitComparison)))
               screened.emplace_back(fit->difference, &candidate);
         }
         std::sort(screened.begin(), screened.end(), [](const auto& a, const auto& b) {return a.first < b.first; });
         std::vector<std::pair<double, const Candidate*>> promising;
         std::set<std::_tstring> promising_factories;
         for (const auto& s : screened)
         {
            if (promising.size() < 5 && promising_factories.insert(s.second->factory->GetName()).second)
               promising.push_back(s);
         }
         count = 0;
         for (const auto& [difference, candidate] : promising)
         {
            std::_tostringstream progress;
            progress << _T("Fitting the section of ") << A2T(beam_name.c_str()) << _T(": full fit ") << ++count << _T(" of ") << promising.size() << _T(" (") << candidate->factory->GetName() << _T(")");
            pProgress->UpdateMessage(progress.str().c_str());
            consider(*candidate, prismatic(candidate->factory, m_pFitter->Fit(candidate->factory, candidate->start, section, *m_pFitComparison)));
         }

         std::ostringstream os_how;
         os_how << "all dimensions fitted, " << promising_factories.size() << " of " << factory_names.size() << " beam factories fully fitted after screening; the closest:";
         for (size_t i = 0; i < promising.size(); i++)
            os_how << (i == 0 ? " " : ", ") << T2A(promising[i].second->factory->GetName().c_str()) << " " << std::fixed << std::setprecision(1) << promising[i].first;
         if (!varying_depth.empty())
         {
            os_how << "; can't make the section as a girder of constant depth:";
            for (const auto& name : varying_depth)
               os_how << " " << name << (&name == &*varying_depth.rbegin() ? "" : ",");
         }
         how = os_how.str();
      }
      auto fit_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - fit_start_time).count();

      // Several beam factories can make the same shape (e.g. the NU girder of the Nebraska NU Girder factory and of other
      // publishers' factories). The fits as good as the best (accepted, of the same rank, with outlines within the length
      // tolerance of each other) from other factories are ties: the user picks the factory, or, without a user (e.g. a
      // command line import), the one published by WSDOT, else the best fit
      if (best_fit && best_fit->differences.empty())
      {
         std::vector<std::pair<const Candidate*, const DimensionFit*>> ties{ { best_candidate, &*best_fit } };
         std::set<std::_tstring> tie_names{ best_candidate->factory->GetName() };
         for (const auto& [candidate, fit] : fits)
         {
            if (fit.differences.empty() && rank(*candidate, fit) == rank(*best_candidate, *best_fit) &&
               tie_names.insert(candidate->factory->GetName()).second &&
               OutlineDistance(fit.section.outline, best_fit->section.outline) <= SectionTolerances().length)
               ties.emplace_back(candidate, &fit);
         }
         // The model's girder has a constant depth: a factory made for such girders (its default girder has the same
         // height at both ends) is preferred to one made for girders that vary in depth (e.g. a tapered I-beam), which
         // only ties by being fitted with a constant depth
         auto made_for_constant_depth = [](const Candidate* candidate) {return HasConstantDepthDefaults(*candidate->factory); };
         if (std::any_of(ties.begin(), ties.end(), [&](const auto& tie) {return made_for_constant_depth(tie.first); }))
         {
            ties.erase(std::remove_if(ties.begin(), ties.end(), [&](const auto& tie) {return !made_for_constant_depth(tie.first); }), ties.end());
            if (std::find_if(ties.begin(), ties.end(), [&](const auto& tie) {return tie.first == best_candidate; }) == ties.end())
            {
               best_candidate = ties.front().first;
               best_fit = *ties.front().second;
            }
         }
         if (1 < ties.size())
         {
            tie_names.clear();
            for (const auto& tie : ties)
               tie_names.insert(tie.first->factory->GetName());

            auto is_wsdot = [](const std::_tstring& publisher)
               {
                  CString p(publisher.c_str());
                  p.MakeUpper();
                  return p == _T("WSDOT") || p.Find(_T("WASHINGTON STATE DEPARTMENT OF TRANSPORTATION")) != -1;
               };
            size_t default_choice = 0;
            for (size_t i = 0; i < ties.size(); i++)
            {
               if (is_wsdot(ties[i].first->factory->GetPublisher()))
               {
                  default_choice = i;
                  break;
               }
            }

            size_t choice = default_choice;
            std::string chosen_by = "published by " + std::string(T2A(ties[default_choice].first->factory->GetPublisher().c_str()));
            auto found = m_FactoryChoices.find(tie_names);
            if (found != m_FactoryChoices.end())
            {
               for (size_t i = 0; i < ties.size(); i++)
               {
                  if (ties[i].first->factory->GetName() == found->second)
                     choice = i;
               }
               chosen_by = "chosen earlier in this import";
            }
            else if (m_Importer.IsInteractive())
            {
               std::_tostringstream responses;
               for (const auto& [candidate, fit] : ties)
               {
                  responses << candidate->factory->GetName() << _T(" (") << candidate->factory->GetGirderFamilyName() << _T(", published by ") << candidate->factory->GetPublisher();
                  if (!candidate->entry.empty())
                     responses << _T(", from ") << candidate->entry;
                  responses << _T(")") << std::endl;
               }
               std::_tostringstream question;
               question << _T("The section of ") << A2T(beam_name.c_str()) << _T(" fits these beam factories equally. Select the beam factory for its library entry.");
               int result = AfxChoose(_T("Select Beam Factory"), question.str().c_str(), responses.str().c_str(), (int)default_choice, FALSE);
               if (0 <= result && result < (int)ties.size())
               {
                  choice = (size_t)result;
                  chosen_by = "chosen by the user";
               }
            }
            m_FactoryChoices[tie_names] = ties[choice].first->factory->GetName();

            os << " These beam factories fit the section equally:";
            for (const auto& [candidate, fit] : ties)
               os << " " << T2A(candidate->factory->GetName().c_str()) << " (" << T2A(candidate->factory->GetPublisher().c_str()) << ")" << (&candidate == &ties.back().first ? ";" : ",");
            os << " " << T2A(ties[choice].first->factory->GetName().c_str()) << " was " << chosen_by << ".";

            best_candidate = ties[choice].first;
            best_fit = *ties[choice].second;
         }
      }

      std::_tstring best_start = best_candidate ? best_candidate->entry : std::_tstring();
      std::string best_factory = best_candidate ? std::string(T2A(best_candidate->factory->GetName().c_str())) : std::string();

      if (best_fit && best_fit->differences.empty())
      {
         std::string grid;
         entry = CreateFittedEntry(best_candidate->factory, girder_type_names, beam_name, *best_fit, grid);
         os << " The " << best_factory << " beam factory (" << T2A(best_candidate->factory->GetGirderFamilyName().c_str()) << ") fits the section ("
            << how << ", from " << (best_start.empty() ? std::string("its default dimensions") : std::string(T2A(best_start.c_str()))) << ", "
            << std::fixed << std::setprecision(1) << fit_seconds << " s), so a new library entry, " << T2A(entry->GetName().c_str()) << ", was created with " << grid << ":";
         for (size_t i = 0; i < best_fit->dimensions.size(); i++)
            os << " " << T2A(best_fit->dimensions[i].first.c_str()) << " " << std::fixed << std::setprecision(3)
               << WBFL::Units::ConvertFromSysUnits(best_fit->dimensions[i].second, WBFL::Units::Measure::Inch) << " in" << (i + 1 < best_fit->dimensions.size() ? "," : ".");
         for (const auto& allowed : best_fit->allowed)
            os << " It differs locally: " << allowed << ".";
         m_LibrarySections[entry->GetName()] = best_fit->section;
      }
      else if (best_fit)
      {
         os << " The best fit, from " << (best_start.empty() ? "the defaults of the " + best_factory + " beam factory" : std::string(T2A(best_start.c_str())) + " (" + best_factory + ")")
            << " (" << how << ", " << std::fixed << std::setprecision(1) << fit_seconds << " s), didn't match the section (";
         for (const auto& d : best_fit->differences)
            os << d << (&d == &best_fit->differences.back() ? "" : ", ");
         os << "). The closest fit:";
         const auto& start_dimensions = best_candidate->start;
         for (size_t i = 0; i < best_fit->dimensions.size(); i++)
         {
            const auto& [name, value] = best_fit->dimensions[i];
            os << " " << T2A(name.c_str()) << " " << std::fixed << std::setprecision(3) << (i < start_dimensions.size() ? WBFL::Units::ConvertFromSysUnits(start_dimensions[i].second, WBFL::Units::Measure::Inch) : 0.0)
               << " -> " << WBFL::Units::ConvertFromSysUnits(value, WBFL::Units::Measure::Inch) << (i + 1 < best_fit->dimensions.size() ? "," : ".");
         }
      }
   }
   WBFL::System::Logger::Info(os.str());
   return entry;
}

bool CIfcBridgeImporter::IsUsableLibraryEntry(const std::_tstring& name)
{
   if (m_CreatedEntries.count(name))
      return true;
   GET_IFACE2(m_Importer.GetBroker(), ILibrary, pLibrary);
   const auto* pEntry = pLibrary->GetGirderEntry(name.c_str());
   return pEntry && !pEntry->IsEditingEnabled();
}

const GirderLibraryEntry* CIfcBridgeImporter::CreateFittedEntry(std::shared_ptr<PGS::Beams::BeamFactory> factory, const std::vector<std::string>& girder_type_names, const std::string& beam_name, const DimensionFit& fit, std::string& grid_description)
{
   USES_CONVERSION;
   GET_IFACE2(m_Importer.GetBroker(), ILibrary, pLibrary);
   auto& library = pLibrary->GetGirderLibrary();

   // the name: the girder type name from the model (else the beam's name) with " (IFC)", so it doesn't conflict with an
   // entry of a master library the project is opened with later (e.g. WSDOT's WF42G for a project imported under another
   // configuration), and " (2)", and so on, if it's taken
   std::_tstring base(A2T((girder_type_names.empty() ? beam_name : girder_type_names.front()).c_str()));
   base += _T(" (IFC)");
   std::_tstring name = base;
   for (int n = 2; pLibrary->GetGirderEntry(name.c_str()); n++)
   {
      std::_tostringstream os;
      os << base << _T(" (") << n << _T(")");
      name = os.str();
   }

   GirderLibraryEntry entry;
   entry.SetBeamFactory(factory); // its default dimensions
   for (const auto& [dimension, value] : fit.dimensions)
      entry.SetDimension(dimension, value, false /*the strand grid is generated below*/);

   auto counts = GenerateStrandGrid(entry, fit.section.outline);
   std::ostringstream os;
   os << "a generated 2 in strand grid (" << counts[0] << " straight, " << counts[1] << " harped, and " << counts[2] << " temporary positions)";
   grid_description = os.str();

   GirderLibraryEntry::GirderEntryDataErrorVec errors;
   entry.ValidateData(&errors);
   for (const auto& error : errors)
      WBFL::System::Logger::Warning("New library entry " + std::string(T2A(name.c_str())) + ": " + std::string(T2A(error.GetErrorMsg().c_str())));

   library.AddEntry(entry, name.c_str());
   library.EnableEditing(name.c_str(), true);
   m_CreatedEntries.insert(name);
   return pLibrary->GetGirderEntry(name.c_str());
}

const GirderLibraryEntry* CIfcBridgeImporter::GetGirderLibraryEntry(IfcSchema::IfcBeam beam, const GirderCenterline& girder)
{
   USES_CONVERSION;

   auto beam_name = beam.Name().value_or("unnamed beam");
   auto girder_type_names = GetGirderTypeNames(beam);

   // the girder's section, at mid-length (end blocks and tapers are at the ends)
   std::optional<SectionMeasures> section;
   if (girder.mesh)
   {
      if (auto outline = CutGirderSection(*girder.mesh, girder.start, girder.end, 0.5))
         section = MeasureSection(*outline);
      else
         WBFL::System::Logger::Info("The section of " + beam_name + " couldn't be cut from its geometry. The girder type is found by name.");
   }

   auto found = std::find_if(m_GirderMatches.begin(), m_GirderMatches.end(), [&](const GirderMatch& match)
      {
         if (match.names != girder_type_names || match.section.has_value() != section.has_value())
            return false;
         return !section || m_pMatchComparison->Compare(*section, *match.section).empty() || SameGirderSection(*section, *match.section);
      });
   if (found != m_GirderMatches.end())
   {
      if (section && found->entry)
         WBFL::System::Logger::Info("Section of " + beam_name + ": the same as an earlier girder's, so " + std::string(CT2A(found->entry->GetName().c_str())) + " is used without fitting again.");
      return found->entry;
   }

   if (girder_type_names.empty())
      CIfcImporter::GetTargetReader().ReportNotFound("girder.type_names", beam, beam_name);

   if (section)
   {
      if (auto entry = MatchLibraryEntry(beam_name, girder_type_names, *section))
      {
         m_GirderMatches.push_back({ girder_type_names, section, entry });
         return entry;
      }
   }

   GET_IFACE2_NOCHECK(m_Importer.GetBroker(), ILibrary, pLibrary); // not used when no name is found
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

   m_GirderMatches.push_back({ girder_type_names, section, girder_library_entry });
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

