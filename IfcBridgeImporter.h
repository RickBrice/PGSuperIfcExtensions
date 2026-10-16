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
#pragma once

#include "IfcImporter.h"
#include "BeamSpacing.h"
#include "GirderSection.h"
#include <IFace\Project.h>
#include <map>
#include <set>

class CIfcBridgeImporter
{
public:
   CIfcBridgeImporter(CIfcImporter& importer);
   CIfcImporter::ImportResult Import(ifcopenshell::file& file,bool bDeriveAlignmentFromDeck);

private:
   CIfcImporter& m_Importer;

   IfcSchema::IfcBridge GetBridge(ifcopenshell::file& file);
   bool IsValidBridge(ifcopenshell::file& file, IfcSchema::IfcBridge bridge);
   bool HasValidGirders(ifcopenshell::file& file, IfcSchema::IfcBridge bridge);
   bool HasValidGirdersByTPF(ifcopenshell::file& file, IfcSchema::IfcBridge bridge);
   bool HasValidGirdersByOther(ifcopenshell::file& file, IfcSchema::IfcBridge bridge);

   void SetGirderProperties(ifcopenshell::file& file, CBridgeDescription2& bridge_desc, const GirderLayout& layout);
   void ImportSlab(ifcopenshell::file& file, CBridgeDescription2& bridge_desc, const GirderLayout& layout, const std::vector<Float64>& pier_stations);
   
   const GirderLibraryEntry* GetGirderLibraryEntry(IfcSchema::IfcBeam beam, const GirderCenterline& girder);
   std::vector<std::string> GetGirderTypeNames(IfcSchema::IfcBeam beam);

   // The library entry whose section matches the girder's section (devdocs/GirderImportPlan.md, G1), or nullptr
   const GirderLibraryEntry* MatchLibraryEntry(const std::string& beam_name, const std::vector<std::string>& girder_type_names, const SectionMeasures& section);

   // A new library entry (in the project's library) cloned from a library entry with fitted dimensions (G1 step 3)
   const GirderLibraryEntry* CreateFittedEntry(const std::_tstring& start_name, const std::vector<std::string>& girder_type_names, const std::string& beam_name, const DimensionFit& fit);

   // the girder types already matched (girder type names, section, library entry or nullptr), so girders of a type are
   // matched and logged once: a girder's section only has to match a previous one within the tolerances
   struct GirderMatch
   {
      std::vector<std::string> names;
      std::optional<SectionMeasures> section;
      const GirderLibraryEntry* entry;
   };
   std::vector<GirderMatch> m_GirderMatches;

   // the sections of the library entries, by name (nullopt if the entry's beam factory can't make one)
   std::map<std::_tstring, std::optional<SectionMeasures>> m_LibrarySections;

   // the beam factory chosen among factories that fit a section equally (e.g. the same NU shape from several), by the
   // names of the factories, so the user is asked once per import
   std::map<std::set<std::_tstring>, std::_tstring> m_FactoryChoices;

   // the strategies (see GirderSection.h): recognizing a library entry, accepting a fitted entry, fitting only the height
   // of an entry whose flanges match, and fitting all of the dimensions (a short fit to screen the beam factories, then
   // a full fit of the most promising)
   std::unique_ptr<ISectionComparison> m_pMatchComparison = std::make_unique<ExactSectionComparison>();
   std::unique_ptr<ISectionComparison> m_pFitComparison = std::make_unique<FittedSectionComparison>();
   std::unique_ptr<IDimensionFitter> m_pHeightFitter = std::make_unique<HeightDimensionFitter>();
   std::unique_ptr<IDimensionFitter> m_pScreenFitter = std::make_unique<OutlineDimensionFitter>(10);
   std::unique_ptr<IDimensionFitter> m_pFitter = std::make_unique<OutlineDimensionFitter>(100);
   bool m_bLibrarySectionsReported = false;

   bool DeriveAlignmentFromDeck(ifcopenshell::file& file);

   void Experiment(ifcopenshell::file& file);
};