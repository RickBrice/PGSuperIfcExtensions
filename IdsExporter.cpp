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

// IdsExporter.cpp : implementation of CIdsExporter
//
// NOTE (build): this translation unit is the only one that includes the generated
// xsd-cxx IDS binding (Schema/ids-binding.hxx). The binding drags in Xerces-C++ DOM
// headers, which clash with the <msxml*.h> ::DOMDocument coclass forward-declared via
// stdafx.h under /permissive-. This file therefore builds with ConformanceMode=false
// (see PGSuperIfcExtensions.vcxproj), matching how F:\ARP\WBFL\Units handles the same
// xsd-cxx 4.2.0 generated code.
//
// What this writes: a project-specific buildingSMART IDS (Information Delivery
// Specification) that asserts every PGSuper-computed design value the IFC exporter
// (IfcExporter.cpp + PropertySets.h / QuantitySets.h / Materials.h /
// USBridge_Classifications.h / Rebar.h) puts into a file for a prestressed concrete
// girder: concrete strengths, jacking force/stress, strand & rebar material
// strengths, section geometry, cambers, quantities, usBridge classifications and
// material associations. Grouped one <specification> per girder segment (plus one
// per distinct material / beam type). All numeric values are emitted in SI base
// units (Pa, m, rad, m2, kg) - PGSuper's system units; an IDS validator normalises
// the IFC value to SI before comparing, so this is correct whether the IFC carries
// SI or US display units.

#include "stdafx.h"
#include "IdsExporter.h"
#include "BeamLabels.h"
#include "IdsBuilder.h"
#include "IfcMappingTable.h"

#include <IFace/Tools.h>
#include <IFace/Project.h>
#include <IFace/Bridge.h>
#include <IFace/Intervals.h>
#include <IFace/PointOfInterest.h>
#include <IFace/AnalysisResults.h>
#include <IFace/DocumentType.h>

#include <EAF/EAFDisplayUnits.h>

#include <PsgLib/BridgeDescription2.h>
#include <psgLib/GirderLabel.h>

#include <LRFD/RebarPool.h>
#include <Materials/PsStrand.h>
#include <Materials/Rebar.h>
#include "SteelSpecifications.h"

#include <MfcTools/Format.h>

#include <Units/Measure.h>
#include <Units/Convert.h>

#include <atlconv.h>
#include <array>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
   using namespace ids_builder;

   // ---- small formatting / conversion helpers ------------------------------------

   // Full-precision plain-decimal text for an <ids:value> number (SI base unit).
   std::string FormatValue(double v)
   {
      std::ostringstream os;
      os << std::defaultfloat << std::setprecision(10) << v;
      return os.str();
   }

   // Concrete strength is meaningful only to 0.1 ksi. Round there, then return the
   // equivalent in system units (Pa) so the emitted requirement matches a model whose
   // author transcribed "4.0 KSI".
   double RoundToTenthKsi_SysUnits(double value_SysUnits)
   {
      double ksi = WBFL::Units::ConvertFromSysUnits(value_SysUnits, WBFL::Units::Measure::KSI);
      ksi = std::floor(ksi * 10.0 + 0.5) / 10.0;
      return WBFL::Units::ConvertToSysUnits(ksi, WBFL::Units::Measure::KSI);
   }

   // e.g. "f'ci = 4.0 KSI (27.58 MPa)"
   std::string StressInstruction(const char* symbol, double value_SysUnits)
   {
      double ksi = WBFL::Units::ConvertFromSysUnits(value_SysUnits, WBFL::Units::Measure::KSI);
      double mpa = WBFL::Units::ConvertFromSysUnits(value_SysUnits, WBFL::Units::Measure::MPa);
      std::ostringstream os;
      os << symbol << " = " << std::fixed << std::setprecision(1) << ksi << " KSI ("
         << std::setprecision(2) << mpa << " MPa)";
      return os.str();
   }

   std::string ForceInstruction(const char* symbol, double value_SysUnits)
   {
      double kip = WBFL::Units::ConvertFromSysUnits(value_SysUnits, WBFL::Units::Measure::Kip);
      double kN = WBFL::Units::ConvertFromSysUnits(value_SysUnits, WBFL::Units::Measure::Newton) / 1000.0;
      std::ostringstream os;
      os << symbol << " = " << std::fixed << std::setprecision(1) << kip << " kip ("
         << std::setprecision(1) << kN << " kN)";
      return os.str();
   }

   std::string LengthInstruction(const char* symbol, double value_SysUnits)
   {
      double ft = WBFL::Units::ConvertFromSysUnits(value_SysUnits, WBFL::Units::Measure::Feet);
      double m = WBFL::Units::ConvertFromSysUnits(value_SysUnits, WBFL::Units::Measure::Meter);
      std::ostringstream os;
      os << symbol << " = " << std::fixed << std::setprecision(3) << ft << " ft ("
         << std::setprecision(4) << m << " m)";
      return os.str();
   }

   std::string SegmentIdentifier(const CSegmentKey& key)
   {
      std::ostringstream os;
      os << "G" << (key.groupIndex + 1) << "-G" << (key.girderIndex + 1) << "-S" << (key.segmentIndex + 1);
      return os.str();
   }

   // The IDS schema restricts <ids:author> to an e-mail-like pattern ([^@]+@[^\.]+\..+).
   bool LooksLikeEmail(const CString& s)
   {
      int at = s.Find(_T('@'));
      if (at <= 0) return false;
      int dot = s.Find(_T('.'), at + 2);
      return 0 < dot && dot < s.GetLength() - 1;
   }

   // ---- IDS tree building -------------------------------------------------------

   IDS::idsValue MinInclusiveValue(const std::string& valueText) { return RestrictionValue("xs:double", "xs:minInclusive", valueText); }

   // A pinned numeric value, honoring Exact vs. Minimum matching.
   IDS::idsValue NumericFacet(const CIdsExportOptions& options, double value_SI)
   {
      if (options.strength_match == CIdsExportOptions::StrengthMatch::Exact)
         return SimpleValue(FormatValue(value_SI));
      return MinInclusiveValue(FormatValue(value_SI));
   }

   // ---- requirement-facet factories --------------------------------------------

   IDS::property NumProperty(const CIdsExportOptions& options, const std::string& pset, const std::string& baseName,
                             const char* ifcDataType, double value_SI, const std::string& instruction = {})
   {
      return PropertyReq(pset, baseName, ifcDataType, NumericFacet(options, value_SI), Card::Required, instruction);
   }

   IDS::property StrProperty(const std::string& pset, const std::string& baseName, const std::string& value)
   {
      return PropertyReq(pset, baseName, "IFCLABEL", SimpleValue(value), Card::Required);
   }

   IDS::property PresenceProperty(const std::string& pset, const std::string& baseName, const char* ifcDataType, Card card = Card::Required)
   {
      return PropertyReq(pset, baseName, ifcDataType, std::nullopt, card);
   }

   // ---- applicability / specification scaffolding ------------------------------

   // ---- gathered design values -------------------------------------------------

   struct StrandTypeValues
   {
      bool        present = false;
      double      pjack = 0.0;           // IfcTendon.TensionForce  (N, total for the group)
      double      jacking_stress = 0.0;  // IfcTendon.PreStress     (Pa)
      bool        any_debonded = false;
      std::string coating;               // "NONE" / "EPOXYCOATED"
      std::string material_name;         // strand IfcMaterial.Name
      double      nominal_diameter = 0.0;
      double      cross_section_area = 0.0;
      double      fy = 0.0;
      double      fpu = 0.0;
      std::string grade_label;           // Pset_MaterialSteel.StructuralGrade
   };

   struct SegmentDesignValues
   {
      CSegmentKey key;
      std::string beam_name;             // IfcBeam.Name (GetIfcBeamName)

      // concrete strengths (Pa)
      double fci = 0.0;   // release  -> ReleaseStrength / FormStrippingStrength / LiftingStrength
      double fc28 = 0.0;  // 28-day   -> TransportationStrength, StrengthClass, Pset_MaterialConcrete
      double fpj = 0.0;   // permanent jacking stress -> InitialTension
      double max_agg_size = 0.0; // m

      // geometry
      double span_length = 0.0; // m
      double slope_angle = 0.0; // rad
      double roll_angle = 0.0;  // rad
      double min_support_length = 0.0; // m (bunk point)

      // cambers (m)
      double initial_camber = 0.0;
      double final_camber = 0.0;
      double screed_camber = 0.0;

      // labels
      std::string type_designation;        // girder family name
      std::string shape_name;              // girder library name
      std::string design_location_number;  // segment label

      std::string concrete_material_name;

      std::array<StrandTypeValues, 3> strands;

      std::string rebar_material_name;     // longitudinal rebar IfcMaterial.Name
   };

   struct StrandMaterial { double fy = 0, fpu = 0, db = 0, area = 0; std::string grade_label; };
   struct RebarMaterial  { double fy = 0, fpu = 0, eu = 0; std::string spec, edition; };
   struct ConcreteMaterial { double fc28 = 0, max_agg_size = 0; };

   SegmentDesignValues GatherSegmentValues(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CSegmentKey& segmentKey,
                                           bool bSpliced, const CIdsExportOptions& options)
   {
      USES_CONVERSION;

      GET_IFACE2(pBroker, IMaterials, pMaterials);
      GET_IFACE2(pBroker, IIntervals, pIntervals);
      GET_IFACE2(pBroker, IStrandGeometry, pStrandGeom);
      GET_IFACE2(pBroker, IBridge, pBridge);
      GET_IFACE2(pBroker, IGirder, pGirder);
      GET_IFACE2(pBroker, IBridgeDescription, pIBridgeDesc);
      GET_IFACE2(pBroker, IPointOfInterest, pPoi);
      GET_IFACE2(pBroker, ICamber, pCamber);
      GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

      SegmentDesignValues sdv;
      sdv.key = segmentKey;
      sdv.beam_name = GetIfcBeamName(segmentKey, bSpliced);

      const IntervalIndexType releaseIntervalIdx = pIntervals->GetPrestressReleaseInterval(segmentKey);
      sdv.fci = pMaterials->GetSegmentFc(segmentKey, releaseIntervalIdx);
      sdv.fc28 = pMaterials->GetSegmentFc28(segmentKey);
      sdv.fpj = pStrandGeom->GetJackingStress(segmentKey, pgsTypes::Permanent);
      sdv.max_agg_size = pMaterials->GetSegmentMaxAggrSize(segmentKey);

      sdv.span_length = pBridge->GetSegmentSpanLength(segmentKey);
      sdv.slope_angle = atan(pBridge->GetSegmentSlope(segmentKey));
      sdv.roll_angle = atan(pGirder->GetOrientation(segmentKey));

      const auto& handling = pIBridgeDesc->GetBridgeDescription()
         ->GetGirderGroup(segmentKey.groupIndex)
         ->GetGirder(segmentKey.girderIndex)
         ->GetSegment(segmentKey.segmentIndex)->HandlingData;
      sdv.min_support_length = std::max(handling.LeadingSupportPoint, handling.TrailingSupportPoint);

      PoiList vPoi;
      pPoi->GetPointsOfInterest(segmentKey, POI_RELEASED_SEGMENT | POI_5L, &vPoi);
      if (!vPoi.empty())
      {
         const pgsPointOfInterest& poiMS = vPoi.front();
         sdv.initial_camber = pCamber->GetInitialCamber(poiMS);
         sdv.final_camber = pCamber->GetExcessCamber(poiMS, pgsTypes::CreepTime::Max);
         sdv.screed_camber = pCamber->GetScreedCamber(poiMS, pgsTypes::CreepTime::Max);
      }

      sdv.type_designation = ToUtf8(pIBridgeDesc->GetBridgeDescription()->GetGirderFamilyName());
      sdv.shape_name = ToUtf8(pIBridgeDesc->GetGirder(segmentKey)->GetGirderName());
      {
         // DesignLocationNumber is always the numeric form ("Span 1, Girder 1"),
         // independent of the ambient alpha/numeric label preference - matches
         // the girder.designation target (IfcTargets.cpp). The RAII guard
         // restores the ambient setting so it doesn't leak into GetIfcBeamName() for
         // later segments.
         pgsAutoGirderLabel autoLabel;
         sdv.design_location_number = ToUtf8(SEGMENT_LABEL(segmentKey));
      }

      {
         std::ostringstream os;
         os << "Precast Concrete, f'c = "
            << T2A((LPCTSTR)(::FormatDimension(sdv.fc28, pDisplayUnits->GetStressUnit())));
         sdv.concrete_material_name = os.str();
      }

      for (int i = 0; i < 3; i++)
      {
         pgsTypes::StrandType strandType = pgsTypes::StrandType(i);
         StrandIndexType n = pStrandGeom->GetStrandCount(segmentKey, strandType);
         if (n == 0) continue;

         const auto* pStrand = pMaterials->GetStrandMaterial(segmentKey, strandType);
         StrandTypeValues& s = sdv.strands[i];
         s.present = true;
         s.pjack = pStrandGeom->GetPjack(segmentKey, strandType);
         s.jacking_stress = pStrandGeom->GetJackingStress(segmentKey, strandType);
         s.coating = (pStrand->GetCoating() == WBFL::Materials::PsStrand::Coating::None) ? "NONE" : "EPOXYCOATED";
         s.material_name = T2A(pStrand->GetName().c_str());
         s.nominal_diameter = pStrand->GetNominalDiameter();
         s.cross_section_area = pStrand->GetNominalArea();
         s.fy = pStrand->GetYieldStrength();
         s.fpu = pStrand->GetUltimateStrength();
         {
            std::ostringstream os;
            os << "ASTM A416 Grade "
               << T2A(WBFL::Materials::PsStrand::GetGrade(pStrand->GetGrade(), true).c_str());
            s.grade_label = os.str();
         }
         for (StrandIndexType strandIdx = 0; strandIdx < n && !s.any_debonded; strandIdx++)
         {
            Float64 dbStart, dbEnd;
            if (pStrandGeom->IsStrandDebonded(segmentKey, strandIdx, strandType, nullptr, &dbStart, &dbEnd))
               s.any_debonded = true;
         }
      }

      {
         WBFL::Materials::Rebar::Type barType;
         WBFL::Materials::Rebar::Grade barGrade;
         pMaterials->GetSegmentLongitudinalRebarMaterial(segmentKey, &barType, &barGrade);
         sdv.rebar_material_name = T2A(WBFL::LRFD::RebarPool::GetMaterialName(barType, barGrade).c_str());
      }

      return sdv;
   }

   // ---- property locations from the mapping table ---------------------------------
   //
   // The IDS checks the design values where the IFC export put them, so the property sets, properties,
   // data types, and classifications come from the same mapping table (devdocs/MappingTablesDesign.md, M4).
   // A requirement for a property the table doesn't export is left out.

   struct Location
   {
      std::string pset;
      std::string name;
      std::string type;   // IFC data type of the value (e.g. IFCPRESSUREMEASURE)
      bool primary;       // the location the importer reads the target from, else the first location
      std::optional<TargetValue> value; // a constant value
   };

   // quantities are checked as the measure they hold
   std::string DataType(const PropertyDeclaration& property)
   {
      if (property.type == "IFCQUANTITYLENGTH") return "IFCLENGTHMEASURE";
      if (property.type == "IFCQUANTITYAREA") return "IFCAREAMEASURE";
      if (property.type == "IFCQUANTITYVOLUME") return "IFCVOLUMEMEASURE";
      if (property.type == "IFCQUANTITYWEIGHT") return "IFCMASSMEASURE";
      if (property.type == "IFCQUANTITYCOUNT") return "IFCCOUNTMEASURE";
      return property.type;
   }

   class CTableLocations
   {
   public:
      CTableLocations(const CIfcMappingTable& table) : m_Table(table) {}

      // where a target is exported for an element role
      std::vector<Location> Target(ElementKind role, PropertyOwner attach, std::string_view target) const
      {
         auto exported = m_Table.GetExportedProperties(role, attach, target);
         auto primary = std::find_if(exported.begin(), exported.end(), [](const auto& e) {return e.property->import; });
         if (primary == exported.end())
            primary = exported.begin();

         std::vector<Location> locations;
         for (auto it = exported.begin(); it != exported.end(); it++)
            locations.push_back({ it->pset->name, it->property->name, DataType(*it->property), it == primary, it->property->value });
         return locations;
      }

      // a property with a constant value, as the table declares it, or nullopt if the table doesn't export it
      std::optional<Location> Constant(ElementKind role, PropertyOwner attach, const char* pset, const char* name) const
      {
         auto exported = m_Table.FindExportedProperty(role, attach, pset, name);
         if (!exported || !exported->property->value)
            return std::nullopt;
         return Location{ exported->pset->name, exported->property->name, DataType(*exported->property), true, exported->property->value };
      }

      // the classification references of an element role: (classification system, identification)
      std::vector<std::pair<std::string, std::string>> Classifications(ElementKind role) const
      {
         std::vector<std::pair<std::string, std::string>> classifications;
         for (const auto* c : m_Table.GetClassifications(role))
            classifications.emplace_back(c->system, c->identification);
         return classifications;
      }

   private:
      const CIfcMappingTable& m_Table;
   };

   std::string ConstantText(const Location& location)
   {
      return std::holds_alternative<std::string>(*location.value) ? std::get<std::string>(*location.value) : FormatTargetValue(*location.value);
   }

   double ConstantNumber(const Location& location)
   {
      const auto& v = *location.value;
      return std::holds_alternative<Float64>(v) ? std::get<Float64>(v) : (std::holds_alternative<Int64>(v) ? (double)std::get<Int64>(v) : 0.0);
   }

   // adds the requirements of the table locations: all of them, or only the primary location or only the others
   enum class Which { All, Primary, Others };
   bool Selected(const Location& location, Which which)
   {
      return which == Which::All || (which == Which::Primary) == location.primary;
   }

   void AddNumbers(IDS::requirements& requirements, const CIdsExportOptions& options, const std::vector<Location>& locations, Which which, double value, const std::string& instruction = {})
   {
      for (const auto& l : locations)
      {
         if (Selected(l, which))
            requirements.property().push_back(NumProperty(options, l.pset, l.name, l.type.c_str(), value, instruction));
      }
   }

   void AddTexts(IDS::requirements& requirements, const std::vector<Location>& locations, const std::string& value)
   {
      for (const auto& l : locations)
         requirements.property().push_back(PropertyReq(l.pset, l.name, l.type.c_str(), SimpleValue(value), Card::Required));
   }

   void AddPresence(IDS::requirements& requirements, const std::vector<Location>& locations, Card card = Card::Required)
   {
      for (const auto& l : locations)
         requirements.property().push_back(PresenceProperty(l.pset, l.name, l.type.c_str(), card));
   }

   void AddConstant(IDS::requirements& requirements, const std::optional<Location>& location)
   {
      if (location)
         requirements.property().push_back(PropertyReq(location->pset, location->name, location->type.c_str(), SimpleValue(ConstantText(*location)), Card::Required));
   }

   void AddClassifications(IDS::requirements& requirements, const CTableLocations& table, ElementKind role)
   {
      for (const auto& [system, identification] : table.Classifications(role))
         requirements.classification().push_back(ClassificationReq(system, identification));
   }

   // ---- per-segment spec builders ---------------------------------------------

   IDS::specificationType BuildBeamSpec(const CIdsExportOptions& options, const CTableLocations& table, CIdsExportOptions::BeamId beamId,
                                        const SegmentDesignValues& sdv, const char* predefinedType)
   {
      IDS::applicabilityType applicability = MakeApplicability("IFCBEAM", predefinedType, "1");

      std::string idAttributeName = "Name";
      std::string idValue = sdv.beam_name;
      if (beamId == CIdsExportOptions::BeamId::GlobalId)
      {
         auto it = options.global_id_by_segment.find(sdv.key);
         if (it != options.global_id_by_segment.end() && !it->second.empty())
         {
            idAttributeName = "GlobalId";
            idValue = it->second;
         }
      }
      AddAttr(applicability, idAttributeName, SimpleValue(idValue));

      IDS::requirements requirements;

      const auto G = ElementKind::Girder;
      const auto O = PropertyOwner::Occurrence;

      // structural checks
      AddClassifications(requirements, table, G);
      requirements.material().push_back(MaterialReq("concrete"));

      // the two design strengths (kept even in minimal mode): where the importer reads them
      const double fciR = RoundToTenthKsi_SysUnits(sdv.fci);
      const double fcR = RoundToTenthKsi_SysUnits(sdv.fc28);
      if (options.include_release_strength || options.include_beam_properties)
         AddNumbers(requirements, options, table.Target(G, O, "girder.fci"), Which::Primary, fciR, StressInstruction("f'ci", fciR));
      if (options.include_transportation_strength || options.include_beam_properties)
         AddNumbers(requirements, options, table.Target(G, O, "girder.fc"), Which::Primary, fcR, StressInstruction("f'c", fcR));

      if (options.include_beam_properties)
      {
         AddConstant(requirements, table.Constant(G, O, "Pset_BeamCommon", "Status"));
         AddNumbers(requirements, options, table.Target(G, O, "girder.span"), Which::All, sdv.span_length, LengthInstruction("Span", sdv.span_length));
         AddNumbers(requirements, options, table.Target(G, O, "girder.slope"), Which::All, sdv.slope_angle);
         AddNumbers(requirements, options, table.Target(G, O, "girder.roll"), Which::All, sdv.roll_angle);

         // the strength class label is display-unit formatted, so presence only
         AddPresence(requirements, table.Target(G, O, "girder.strength_class"));

         // the other locations of f'ci (e.g. form stripping and lifting strengths)
         AddNumbers(requirements, options, table.Target(G, O, "girder.fci"), Which::Others, fciR, StressInstruction("f'ci", fciR));
         AddNumbers(requirements, options, table.Target(G, O, "girder.jacking_stress"), Which::All, sdv.fpj, StressInstruction("fpj", sdv.fpj));
         AddTexts(requirements, table.Target(G, O, "girder.family_name"), sdv.type_designation);
         AddNumbers(requirements, options, table.Target(G, O, "girder.bunk_point"), Which::All, sdv.min_support_length, LengthInstruction("bunk point", sdv.min_support_length));
         AddTexts(requirements, table.Target(G, O, "girder.designation"), sdv.design_location_number);

         // populated only when the IFC was exported with camber / batter -> optional
         AddPresence(requirements, table.Target(G, O, "girder.camber_ratio"), Card::Optional);
         AddPresence(requirements, table.Target(G, O, "girder.batter"), Card::Optional);

         AddNumbers(requirements, options, table.Target(G, O, "girder.camber_at_release"), Which::All, sdv.initial_camber, LengthInstruction("D at release", sdv.initial_camber));
         AddNumbers(requirements, options, table.Target(G, O, "girder.camber_after_losses"), Which::All, sdv.final_camber, LengthInstruction("excess camber", sdv.final_camber));
         AddNumbers(requirements, options, table.Target(G, O, "girder.screed_camber"), Which::All, sdv.screed_camber, LengthInstruction("screed camber", sdv.screed_camber));
         AddTexts(requirements, table.Target(G, O, "girder.type_names"), sdv.shape_name);
      }

      if (options.include_quantities)
      {
         // IfcQuantity* values are not reliably unit-normalised by validators -> presence only.
         for (const auto* target : { "girder.length", "girder.cross_section_area", "girder.outer_surface_area", "girder.gross_weight" })
            AddPresence(requirements, table.Target(G, O, target));
      }

      std::string name = ToUtf8(options.specification_name_prefix);
      if (!name.empty()) name += " - ";
      name += sdv.beam_name;

      return MakeSpec(std::move(applicability), name, SegmentIdentifier(sdv.key) + "-BEAM",
         "Precast concrete girder design values from PGSuper", std::move(requirements));
   }

   // ---- bridge-wide (not per-segment) structural/classification specs --------
   //
   // Classification, material and PredefinedType are identical for every occurrence of
   // a kind regardless of which segment it belongs to, so these apply to every matching
   // occurrence in the whole model - entity(+predefinedType) alone, no Name facet needed.

   IDS::specificationType BuildGlobalTendonSpec(const CTableLocations& table)
   {
      IDS::applicabilityType applicability = MakeApplicability("IFCTENDON", "STRAND", "1");

      IDS::requirements requirements;
      AddClassifications(requirements, table, ElementKind::Tendon);
      requirements.material().push_back(MaterialReq("steel"));
      // Exact jacking force/stress values are per (segment, strand type) and can't be
      // asserted without an unambiguous per-segment identifier on the IfcTendon
      // occurrence - just require the attributes are populated.
      requirements.attribute().push_back(AttributeReq("TensionForce", std::nullopt));
      requirements.attribute().push_back(AttributeReq("PreStress", std::nullopt));

      return MakeSpec(std::move(applicability), "Tendon classification and material", "TENDON",
         "Every prestressing strand is a classified steel tendon with jacking data", std::move(requirements));
   }

   IDS::specificationType BuildGlobalTendonBundleSpec(const CTableLocations& table)
   {
      IDS::applicabilityType applicability = MakeApplicability("IFCELEMENTASSEMBLY", "USERDEFINED", "1");

      IDS::requirements requirements;
      AddClassifications(requirements, table, ElementKind::TendonBundle);

      return MakeSpec(std::move(applicability), "Tendon bundle classification", "STRANDS",
         "Every strand bundle is a classified tendon bundle", std::move(requirements));
   }

   IDS::specificationType BuildGlobalRebarSpec(const CTableLocations& table)
   {
      IDS::applicabilityType applicability = MakeApplicability("IFCREINFORCINGBAR", nullptr, "0");

      IDS::requirements requirements;
      AddClassifications(requirements, table, ElementKind::Rebar);
      requirements.material().push_back(MaterialReq("steel"));

      return MakeSpec(std::move(applicability), "Reinforcing bar classification and material", "REBAR",
         "Every reinforcing bar is a classified steel bar", std::move(requirements));
   }

   IDS::specificationType BuildGlobalCageSpec(const CTableLocations& table)
   {
      IDS::applicabilityType applicability = MakeApplicability("IFCELEMENTASSEMBLY", "REINFORCEMENT_UNIT", "0");

      IDS::requirements requirements;
      AddClassifications(requirements, table, ElementKind::ReinforcementCage);

      return MakeSpec(std::move(applicability), "Reinforcement cage classification", "CAGE",
         "Every reinforcement cage is classified", std::move(requirements));
   }

   // ---- bridge-level (per distinct material / beam type) spec builders --------

   IDS::specificationType BuildConcreteMaterialSpec(const CIdsExportOptions& options, const CTableLocations& table, const std::string& name,
                                                    const ConcreteMaterial& m)
   {
      IDS::applicabilityType applicability = MakeApplicability("IFCMATERIAL", nullptr, "1");
      AddAttr(applicability, "Name", SimpleValue(name));

      IDS::requirements requirements;
      const auto G = ElementKind::Girder;
      const auto M = PropertyOwner::Material;
      const double fcR = RoundToTenthKsi_SysUnits(m.fc28);
      AddNumbers(requirements, options, table.Target(G, M, "girder.fc"), Which::All, fcR, StressInstruction("f'c", fcR));
      AddNumbers(requirements, options, table.Target(G, M, "girder.max_aggregate_size"), Which::All, m.max_agg_size, LengthInstruction("max aggregate", m.max_agg_size));

      return MakeSpec(std::move(applicability), "Concrete material - " + name, {},
         "Concrete compressive strength from PGSuper", std::move(requirements));
   }

   IDS::specificationType BuildStrandMaterialSpec(const CIdsExportOptions& options, const CTableLocations& table, const std::string& name,
                                                  const StrandMaterial& m)
   {
      IDS::applicabilityType applicability = MakeApplicability("IFCMATERIAL", nullptr, "1");
      AddAttr(applicability, "Name", SimpleValue(name));

      IDS::requirements requirements;
      const auto T = ElementKind::Tendon;
      const auto M = PropertyOwner::Material;
      AddNumbers(requirements, options, table.Target(T, M, "tendon.fy"), Which::All, m.fy, StressInstruction("fpy", m.fy));
      AddNumbers(requirements, options, table.Target(T, M, "tendon.fpu"), Which::All, m.fpu, StressInstruction("fpu", m.fpu));
      if (auto strain = table.Constant(T, M, "Pset_MaterialSteel", "UltimateStrain"))
         requirements.property().push_back(NumProperty(options, strain->pset, strain->name, strain->type.c_str(), ConstantNumber(*strain)));
      AddTexts(requirements, table.Target(T, M, "tendon.grade"), m.grade_label);
      AddConstant(requirements, table.Constant(T, M, "usBrPset_ACI_TendonMaterial", "Specification"));

      return MakeSpec(std::move(applicability), "Strand material - " + name, {},
         "Prestressing strand material properties from PGSuper", std::move(requirements));
   }

   IDS::specificationType BuildRebarMaterialSpec(const CIdsExportOptions& options, const CTableLocations& table, const std::string& name,
                                                 const RebarMaterial& m)
   {
      IDS::applicabilityType applicability = MakeApplicability("IFCMATERIAL", nullptr, "0");
      AddAttr(applicability, "Name", SimpleValue(name));

      IDS::requirements requirements;
      const auto R = ElementKind::Rebar;
      const auto M = PropertyOwner::Material;
      AddNumbers(requirements, options, table.Target(R, M, "rebar.fy"), Which::All, m.fy, StressInstruction("fy", m.fy));
      AddNumbers(requirements, options, table.Target(R, M, "rebar.fu"), Which::All, m.fpu, StressInstruction("fu", m.fpu));
      AddNumbers(requirements, options, table.Target(R, M, "rebar.elongation"), Which::All, m.eu);
      AddTexts(requirements, table.Target(R, M, "rebar.grade"), name);
      AddTexts(requirements, table.Target(R, M, "rebar.specification"), m.spec);
      AddTexts(requirements, table.Target(R, M, "rebar.specification_edition"), m.edition);

      return MakeSpec(std::move(applicability), "Reinforcing material - " + name, {},
         "Reinforcing steel material properties from PGSuper", std::move(requirements));
   }

   IDS::specificationType BuildTendonTypeSpec(const CIdsExportOptions& options, const std::string& name,
                                              const StrandMaterial& m)
   {
      IDS::applicabilityType applicability = MakeApplicability("IFCTENDONTYPE", "STRAND", "1");
      AddAttr(applicability, "Name", SimpleValue(name));

      IDS::requirements requirements;
      requirements.attribute().push_back(AttributeReq("NominalDiameter", NumericFacet(options, m.db),
         LengthInstruction("strand diameter", m.db)));
      requirements.attribute().push_back(AttributeReq("CrossSectionArea", NumericFacet(options, m.area)));

      return MakeSpec(std::move(applicability), "Strand type - " + name, {},
         "Nominal strand diameter and area from PGSuper", std::move(requirements));
   }

   IDS::specificationType BuildBeamTypeSpec(const CTableLocations& table, const std::string& name)
   {
      IDS::applicabilityType applicability = MakeApplicability("IFCBEAMTYPE", nullptr, "1");
      AddAttr(applicability, "Name", SimpleValue(name));

      IDS::requirements requirements;
      AddConstant(requirements, table.Constant(ElementKind::Girder, PropertyOwner::Type, "Pset_ConcreteElementGeneral", "AssemblyPlace"));
      AddConstant(requirements, table.Constant(ElementKind::Girder, PropertyOwner::Type, "Pset_ConcreteElementGeneral", "CastingMethod"));

      return MakeSpec(std::move(applicability), "Precast girder type - " + name, {},
         "Girder is factory-precast concrete", std::move(requirements));
   }

} // anonymous namespace

bool CIdsExporter::BuildSpecification(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIdsExportOptions& options, const CString& strFilePath)
{
   USES_CONVERSION;
   std::ofstream ofs(T2A(strFilePath), std::ios::binary);
   if (!ofs.good()) return false;
   return BuildSpecification(pBroker, options, ofs);
}

bool CIdsExporter::BuildSpecification(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIdsExportOptions& options, std::ostream& os)
{
   USES_CONVERSION;

   try
   {
      XercesGuard xercesGuard;

      // the property locations of the IFC export. Throws CIfcMappingTableException if the table can't be used
      std::filesystem::path table_path(options.mapping_file.GetString());
      auto pTable = CIfcMappingTable::LoadActive(table_path);
      CTableLocations table(*pTable);

      GET_IFACE2(pBroker, IDocumentType, pDocType);
      GET_IFACE2(pBroker, IBridgeDescription, pIBridgeDesc);
      GET_IFACE2(pBroker, IMaterials, pMaterials);
      GET_IFACE2_NOCHECK(pBroker, IProjectProperties, pProjectProperties);

      const bool bSpliced = pDocType->IsPGSpliceDocument();
      const char* const predefinedType = bSpliced ? "GIRDER_SEGMENT" : "BEAM";

      // GlobalId identification is only usable when the caller has supplied the beam
      // GlobalIds (PGSuperDataExporter.cpp always does); fall back to Name otherwise.
      const CIdsExportOptions::BeamId effectiveBeamId =
         (options.beam_id == CIdsExportOptions::BeamId::GlobalId && options.global_id_by_segment.empty())
            ? CIdsExportOptions::BeamId::Name
            : options.beam_id;

      // ---- <ids:info> ----
      CString title = options.title;
      if (title.IsEmpty()) title = pProjectProperties->GetBridgeName();
      if (title.IsEmpty()) title = _T("PGSuper Bridge");

      IDS::info info(ToUtf8(title));
      if (!options.copyright.IsEmpty())   info.copyright(ToUtf8(options.copyright));
      if (!options.version.IsEmpty())     info.version(ToUtf8(options.version));

      CString description = options.description;
      if (description.IsEmpty())
      {
         description.Format(_T("Prestressed concrete girder design values from PGSuper (%s)."),
                            bSpliced ? _T("PGSplice document") : _T("PGSuper document"));
      }
      info.description(ToUtf8(description));

      if (LooksLikeEmail(options.author)) info.author(IDS::author(ToUtf8(options.author)));
      info.date(Today());
      if (!options.purpose.IsEmpty())     info.purpose(ToUtf8(options.purpose));
      if (!options.milestone.IsEmpty())   info.milestone(ToUtf8(options.milestone));

      IDS::specificationsType specifications;

      // distinct materials / beam types, emitted once after the segment loop
      std::map<std::string, ConcreteMaterial> concreteMaterials;
      std::map<std::string, StrandMaterial>   strandMaterials; // also drives strand-type specs
      std::map<std::string, RebarMaterial>    rebarMaterials;
      std::set<std::string>                   beamTypeNames;
      bool anyStrandsAnywhere = false;
      bool anyRebarAnywhere = false;

      const GroupIndexType nGroups = pIBridgeDesc->GetGirderGroupCount();
      for (GroupIndexType grpIdx = 0; grpIdx < nGroups; grpIdx++)
      {
         const GirderIndexType nGirders = pIBridgeDesc->GetGirderGroup(grpIdx)->GetGirderCount();
         for (GirderIndexType gdrIdx = 0; gdrIdx < nGirders; gdrIdx++)
         {
            const SegmentIndexType nSegments = pIBridgeDesc->GetGirderGroup(grpIdx)->GetGirder(gdrIdx)->GetSegmentCount();
            for (SegmentIndexType segIdx = 0; segIdx < nSegments; segIdx++)
            {
               const CSegmentKey segmentKey(grpIdx, gdrIdx, segIdx);
               const SegmentDesignValues sdv = GatherSegmentValues(pBroker, segmentKey, bSpliced, options);

               beamTypeNames.insert(bSpliced ? std::string("Precast Girder Type") : sdv.shape_name);
               concreteMaterials[sdv.concrete_material_name] = ConcreteMaterial{ sdv.fc28, sdv.max_agg_size };

               specifications.specification().push_back(
                  BuildBeamSpec(options, table, effectiveBeamId, sdv, predefinedType));

               if (options.include_strands)
               {
                  for (int i = 0; i < 3; i++)
                  {
                     const StrandTypeValues& s = sdv.strands[i];
                     if (!s.present) continue;
                     anyStrandsAnywhere = true;
                     strandMaterials[s.material_name] = StrandMaterial{
                        s.fy, s.fpu, s.nominal_diameter, s.cross_section_area, s.grade_label };
                  }
               }

               if (options.include_rebar)
               {
                  // fy/fpu/eu are bar-size independent; use #3 like the IFC exporter.
                  WBFL::Materials::Rebar::Type barType;
                  WBFL::Materials::Rebar::Grade barGrade;
                  pMaterials->GetSegmentLongitudinalRebarMaterial(segmentKey, &barType, &barGrade);
                  const auto* pRebar = WBFL::LRFD::RebarPool::GetInstance()->GetRebar(
                     barType, barGrade, WBFL::Materials::Rebar::Size::bs3);
                  if (pRebar != nullptr)
                  {
                     anyRebarAnywhere = true;
                     rebarMaterials[sdv.rebar_material_name] = RebarMaterial{
                        pRebar->GetYieldStrength(), pRebar->GetUltimateStrength(), pRebar->GetElongation(),
                        GetRebarSpecification(pRebar), GetRebarSpecificationEdition(pRebar) };
                  }
               }
            }
         }
      }

      if (options.include_beam_properties)
      {
         for (const auto& name : beamTypeNames)
            specifications.specification().push_back(BuildBeamTypeSpec(table, name));
      }

      if (options.include_concrete_material)
      {
         for (const auto& [name, m] : concreteMaterials)
            specifications.specification().push_back(BuildConcreteMaterialSpec(options, table, name, m));
      }

      if (options.include_strands)
      {
         for (const auto& [name, m] : strandMaterials)
         {
            specifications.specification().push_back(BuildStrandMaterialSpec(options, table, name, m));
            specifications.specification().push_back(BuildTendonTypeSpec(options, name, m));
         }
         if (anyStrandsAnywhere)
         {
            specifications.specification().push_back(BuildGlobalTendonSpec(table));
            specifications.specification().push_back(BuildGlobalTendonBundleSpec(table));
         }
      }

      if (options.include_rebar)
      {
         for (const auto& [name, m] : rebarMaterials)
            specifications.specification().push_back(BuildRebarMaterialSpec(options, table, name, m));
         if (anyRebarAnywhere)
         {
            specifications.specification().push_back(BuildGlobalRebarSpec(table));
            specifications.specification().push_back(BuildGlobalCageSpec(table));
         }
      }

      if (specifications.specification().empty()) return false;

      IDS::ids document(info, specifications);

      xml_schema::namespace_infomap map;
      map["ids"].name = IDS_NAMESPACE;
      map["ids"].schema = IDS_SCHEMA_LOCATION;
      map["xs"].name = XS_NAMESPACE;

      IDS::ids_(os, document, map, "UTF-8");
      return os.good();
   }
   catch (const CIfcMappingTableException&)
   {
      throw; // the caller reports the message (e.g. what to do about a missing table)
   }
   catch (const xml_schema::exception&)
   {
      return false;
   }
   catch (...)
   {
      return false;
   }
}
