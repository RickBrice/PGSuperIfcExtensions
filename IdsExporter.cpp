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

#include "Schema/ids-binding.hxx"

#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>

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
   const char* const IDS_NAMESPACE = "http://standards.buildingsmart.org/IDS";
   const char* const IDS_SCHEMA_LOCATION = "http://standards.buildingsmart.org/IDS/1.0/ids.xsd";
   const char* const XS_NAMESPACE = "http://www.w3.org/2001/XMLSchema";

   // ---- small formatting / conversion helpers ------------------------------------

   std::string ToUtf8(LPCTSTR s)
   {
      if (s == nullptr) return std::string();
      CW2A conv(s, CP_UTF8);
      return std::string((LPCSTR)conv);
   }

   std::string ToUtf8(const CString& s) { return ToUtf8((LPCTSTR)s); }

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

   xml_schema::date Today()
   {
      std::time_t t = std::time(nullptr);
      std::tm lt{};
      localtime_s(&lt, &t);
      return xml_schema::date(lt.tm_year + 1900,
                              static_cast<unsigned short>(lt.tm_mon + 1),
                              static_cast<unsigned short>(lt.tm_mday));
   }

   // The IDS schema restricts <ids:author> to an e-mail-like pattern ([^@]+@[^\.]+\..+).
   bool LooksLikeEmail(const CString& s)
   {
      int at = s.Find(_T('@'));
      if (at <= 0) return false;
      int dot = s.Find(_T('.'), at + 2);
      return 0 < dot && dot < s.GetLength() - 1;
   }

   // Reproduce (locally - PropertySets.h is IFC-schema-templated and can't be included
   // here) the specification strings the IFC exporter writes.
   std::string RebarSpecification(const WBFL::Materials::Rebar* pRebar)
   {
      switch (pRebar->GetType())
      {
      case WBFL::Materials::Rebar::Type::A615:  return "ASTM A615 (AASHTO M31)";
      case WBFL::Materials::Rebar::Type::A706:  return "ASTM A706";
      case WBFL::Materials::Rebar::Type::A1035: return "ASTM A1035";
      default:                                  return "Unknown";
      }
   }

   std::string RebarSpecificationEdition(const WBFL::Materials::Rebar* pRebar)
   {
      switch (pRebar->GetType())
      {
      case WBFL::Materials::Rebar::Type::A615:  return "2026";
      case WBFL::Materials::Rebar::Type::A706:  return "2026";
      case WBFL::Materials::Rebar::Type::A1035: return "2024";
      default:                                  return "Unknown";
      }
   }

   // ---- IDS tree building -------------------------------------------------------

   IDS::idsValue SimpleValue(const std::string& text)
   {
      IDS::idsValue v;
      v.simpleValue(text);
      return v;
   }

   // Xerces-C++ requires initialization before any DOM/serialization use. Initialize()
   // is reference counted, so pairing it with Terminate() is safe whether or not the
   // host has already initialized the library.
   struct XercesGuard
   {
      XercesGuard() { xercesc::XMLPlatformUtils::Initialize(); }
      ~XercesGuard() { xercesc::XMLPlatformUtils::Terminate(); }
      XercesGuard(const XercesGuard&) = delete;
      XercesGuard& operator=(const XercesGuard&) = delete;
   };

   // RAII helper for a transcoded XMLCh* string.
   struct XStr
   {
      XMLCh* p;
      explicit XStr(const char* s) : p(xercesc::XMLString::transcode(s)) {}
      ~XStr() { xercesc::XMLString::release(&p); }
      operator const XMLCh* () const { return p; }
      XStr(const XStr&) = delete;
      XStr& operator=(const XStr&) = delete;
   };

   // <ids:value><xs:restriction base="..."><xs:<facet> value="..."/></xs:restriction></ids:value>
   IDS::idsValue RestrictionValue(const char* base, const char* facet, const std::string& facetValue)
   {
      IDS::idsValue v;
      xercesc::DOMDocument& doc = v.dom_document();

      XStr xsNs(XS_NAMESPACE);

      xercesc::DOMElement* restriction = doc.createElementNS(xsNs, XStr("xs:restriction"));
      restriction->setAttribute(XStr("base"), XStr(base));

      xercesc::DOMElement* facetElem = doc.createElementNS(xsNs, XStr(facet));
      facetElem->setAttribute(XStr("value"), XStr(facetValue.c_str()));
      restriction->appendChild(facetElem);

      v.any(restriction);
      return v;
   }

   IDS::idsValue MinInclusiveValue(const std::string& valueText) { return RestrictionValue("xs:double", "xs:minInclusive", valueText); }

   // A pinned numeric value, honoring Exact vs. Minimum matching.
   IDS::idsValue NumericFacet(const CIdsExportOptions& options, double value_SI)
   {
      if (options.strength_match == CIdsExportOptions::StrengthMatch::Exact)
         return SimpleValue(FormatValue(value_SI));
      return MinInclusiveValue(FormatValue(value_SI));
   }

   enum class Card { Required, Optional };

   const char* CardText(Card c) { return c == Card::Optional ? "optional" : "required"; }

   // ---- requirement-facet factories --------------------------------------------

   IDS::property PropertyReq(const std::string& pset, const std::string& baseName, const char* ifcDataType,
                             std::optional<IDS::idsValue> value, Card card, const std::string& instruction = {})
   {
      IDS::property p(SimpleValue(pset), SimpleValue(baseName));
      if (ifcDataType) p.dataType(IDS::upperCaseName(ifcDataType));
      p.cardinality(IDS::conditionalCardinality(CardText(card)));
      if (!instruction.empty()) p.instructions(instruction);
      if (value) p.value(std::move(*value));
      return p;
   }

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

   IDS::classification ClassificationReq(const std::string& code)
   {
      IDS::classification c(SimpleValue("usBridge"));
      c.value(SimpleValue(code));           // matches IfcClassificationReference.Identification
      c.cardinality(IDS::conditionalCardinality("required"));
      return c;
   }

   IDS::material MaterialReq(const std::string& nameOrCategory)
   {
      IDS::material m;
      m.value(SimpleValue(nameOrCategory)); // matched against IfcMaterial.Name and .Category
      m.cardinality(IDS::conditionalCardinality("required"));
      return m;
   }

   IDS::attribute AttributeReq(const std::string& name, std::optional<IDS::idsValue> value, const std::string& instruction = {})
   {
      IDS::attribute a(SimpleValue(name));
      a.cardinality(IDS::conditionalCardinality("required"));
      if (!instruction.empty()) a.instructions(instruction);
      if (value) a.value(std::move(*value));
      return a;
   }

   // ---- applicability / specification scaffolding ------------------------------

   IDS::applicabilityType MakeApplicability(const char* ifcClass, const char* predefinedType, const std::string& minOccurs)
   {
      IDS::entityType entity(SimpleValue(ifcClass));
      if (predefinedType) entity.predefinedType(SimpleValue(predefinedType));

      IDS::applicabilityType applicability;
      applicability.entity(entity);
      applicability.minOccurs(minOccurs);
      applicability.maxOccurs(std::string("unbounded"));
      return applicability;
   }

   void AddAttr(IDS::applicabilityType& applicability, const std::string& name, IDS::idsValue value)
   {
      IDS::attributeType attribute(SimpleValue(name));
      attribute.value(std::move(value));
      applicability.attribute().push_back(attribute);
   }

   IDS::specificationType MakeSpec(IDS::applicabilityType applicability, const std::string& name,
                                   const std::string& identifier, const std::string& description,
                                   IDS::requirements requirements)
   {
      IDS::ifcVersion ifcVersion;
      ifcVersion.push_back(IDS::ifcVersion_item("IFC4X3_ADD2"));

      IDS::specificationType spec(applicability, name, ifcVersion);
      if (!identifier.empty())  spec.identifier(identifier);
      if (!description.empty()) spec.description(description);
      spec.requirements(std::move(requirements));
      return spec;
   }

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
         // Create_Pset_PrecastConcreteElementGeneral (PropertySets.h). The RAII guard
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

   // ---- per-segment spec builders ---------------------------------------------

   IDS::specificationType BuildBeamSpec(const CIdsExportOptions& options, CIdsExportOptions::BeamId beamId,
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

      // structural checks
      requirements.classification().push_back(ClassificationReq("usBridge_GirderPrecastConcrete"));
      requirements.material().push_back(MaterialReq("concrete"));

      // the two design strengths (kept even in minimal mode)
      if (options.include_release_strength || options.include_beam_properties)
         requirements.property().push_back(NumProperty(options, "Pset_PrecastConcreteElementGeneral", "ReleaseStrength",
            "IFCPRESSUREMEASURE", RoundToTenthKsi_SysUnits(sdv.fci), StressInstruction("f'ci", RoundToTenthKsi_SysUnits(sdv.fci))));
      if (options.include_transportation_strength || options.include_beam_properties)
         requirements.property().push_back(NumProperty(options, "Pset_PrecastConcreteElementGeneral", "TransportationStrength",
            "IFCPRESSUREMEASURE", RoundToTenthKsi_SysUnits(sdv.fc28), StressInstruction("f'c", RoundToTenthKsi_SysUnits(sdv.fc28))));

      if (options.include_beam_properties)
      {
         const double fciR = RoundToTenthKsi_SysUnits(sdv.fci);

         // Pset_BeamCommon
         requirements.property().push_back(PropertyReq("Pset_BeamCommon", "Status", "IFCLABEL", SimpleValue("NEW"), Card::Required));
         requirements.property().push_back(NumProperty(options, "Pset_BeamCommon", "Span", "IFCPOSITIVELENGTHMEASURE",
            sdv.span_length, LengthInstruction("Span", sdv.span_length)));
         requirements.property().push_back(NumProperty(options, "Pset_BeamCommon", "Slope", "IFCPLANEANGLEMEASURE", sdv.slope_angle));
         requirements.property().push_back(NumProperty(options, "Pset_BeamCommon", "Roll", "IFCPLANEANGLEMEASURE", sdv.roll_angle));

         // Pset_ConcreteElementGeneral (beam-level override) - the label is display-unit
         // formatted, so presence only.
         requirements.property().push_back(PresenceProperty("Pset_ConcreteElementGeneral", "StrengthClass", "IFCLABEL"));

         // Pset_PrecastConcreteElementGeneral
         requirements.property().push_back(NumProperty(options, "Pset_PrecastConcreteElementGeneral", "FormStrippingStrength",
            "IFCPRESSUREMEASURE", fciR, StressInstruction("f'ci", fciR)));
         requirements.property().push_back(NumProperty(options, "Pset_PrecastConcreteElementGeneral", "LiftingStrength",
            "IFCPRESSUREMEASURE", fciR, StressInstruction("f'ci", fciR)));
         requirements.property().push_back(NumProperty(options, "Pset_PrecastConcreteElementGeneral", "InitialTension",
            "IFCPRESSUREMEASURE", sdv.fpj, StressInstruction("fpj", sdv.fpj)));
         requirements.property().push_back(StrProperty("Pset_PrecastConcreteElementGeneral", "TypeDesignation", sdv.type_designation));
         requirements.property().push_back(NumProperty(options, "Pset_PrecastConcreteElementGeneral", "MinimumAllowableSupportLength",
            "IFCPOSITIVELENGTHMEASURE", sdv.min_support_length, LengthInstruction("bunk point", sdv.min_support_length)));
         requirements.property().push_back(StrProperty("Pset_PrecastConcreteElementGeneral", "DesignLocationNumber", sdv.design_location_number));
         // populated only when the IFC was exported with camber / batter -> optional
         requirements.property().push_back(PresenceProperty("Pset_PrecastConcreteElementGeneral", "CamberAtMidspan", "IFCRATIOMEASURE", Card::Optional));
         requirements.property().push_back(PresenceProperty("Pset_PrecastConcreteElementGeneral", "BatterAtStart", "IFCPLANEANGLEMEASURE", Card::Optional));
         requirements.property().push_back(PresenceProperty("Pset_PrecastConcreteElementGeneral", "BatterAtEnd", "IFCPLANEANGLEMEASURE", Card::Optional));

         // usBrPset_PrecastConcreteBeam
         requirements.property().push_back(NumProperty(options, "usBrPset_PrecastConcreteBeam", "MidSpanCamberAtRelease",
            "IFCLENGTHMEASURE", sdv.initial_camber, LengthInstruction("D at release", sdv.initial_camber)));
         requirements.property().push_back(NumProperty(options, "usBrPset_PrecastConcreteBeam", "MidSpanCamberAfterLosses",
            "IFCLENGTHMEASURE", sdv.final_camber, LengthInstruction("excess camber", sdv.final_camber)));
         requirements.property().push_back(NumProperty(options, "usBrPset_PrecastConcreteBeam", "DeflectionShortTerm",
            "IFCLENGTHMEASURE", sdv.screed_camber, LengthInstruction("screed camber", sdv.screed_camber)));
         requirements.property().push_back(StrProperty("usBrPset_PrecastConcreteBeam", "ShapeName", sdv.shape_name));
      }

      if (options.include_quantities)
      {
         // IfcQuantity* values are not reliably unit-normalised by validators -> presence only.
         requirements.property().push_back(PresenceProperty("Qto_BeamBaseQuantities", "Length", "IFCLENGTHMEASURE"));
         requirements.property().push_back(PresenceProperty("Qto_BeamBaseQuantities", "CrossSectionArea", "IFCAREAMEASURE"));
         requirements.property().push_back(PresenceProperty("Qto_BeamBaseQuantities", "OuterSurfaceArea", "IFCAREAMEASURE"));
         requirements.property().push_back(PresenceProperty("Qto_BeamBaseQuantities", "GrossWeight", "IFCMASSMEASURE"));
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

   IDS::specificationType BuildGlobalTendonSpec()
   {
      IDS::applicabilityType applicability = MakeApplicability("IFCTENDON", "STRAND", "1");

      IDS::requirements requirements;
      requirements.classification().push_back(ClassificationReq("usBridge_Tendon"));
      requirements.material().push_back(MaterialReq("steel"));
      // Exact jacking force/stress values are per (segment, strand type) and can't be
      // asserted without an unambiguous per-segment identifier on the IfcTendon
      // occurrence - just require the attributes are populated.
      requirements.attribute().push_back(AttributeReq("TensionForce", std::nullopt));
      requirements.attribute().push_back(AttributeReq("PreStress", std::nullopt));

      return MakeSpec(std::move(applicability), "Tendon classification and material", "TENDON",
         "Every prestressing strand is a classified steel tendon with jacking data", std::move(requirements));
   }

   IDS::specificationType BuildGlobalTendonBundleSpec()
   {
      IDS::applicabilityType applicability = MakeApplicability("IFCELEMENTASSEMBLY", "USERDEFINED", "1");

      IDS::requirements requirements;
      requirements.classification().push_back(ClassificationReq("usBridge_TendonBundle"));

      return MakeSpec(std::move(applicability), "Tendon bundle classification", "STRANDS",
         "Every strand bundle is a classified tendon bundle", std::move(requirements));
   }

   IDS::specificationType BuildGlobalRebarSpec()
   {
      IDS::applicabilityType applicability = MakeApplicability("IFCREINFORCINGBAR", nullptr, "0");

      IDS::requirements requirements;
      requirements.classification().push_back(ClassificationReq("usBridge_ReinforcingBar"));
      requirements.material().push_back(MaterialReq("steel"));

      return MakeSpec(std::move(applicability), "Reinforcing bar classification and material", "REBAR",
         "Every reinforcing bar is a classified steel bar", std::move(requirements));
   }

   IDS::specificationType BuildGlobalCageSpec()
   {
      IDS::applicabilityType applicability = MakeApplicability("IFCELEMENTASSEMBLY", "REINFORCEMENT_UNIT", "0");

      IDS::requirements requirements;
      requirements.classification().push_back(ClassificationReq("usBridge_ReinforcementCage"));

      return MakeSpec(std::move(applicability), "Reinforcement cage classification", "CAGE",
         "Every reinforcement cage is classified", std::move(requirements));
   }

   // ---- bridge-level (per distinct material / beam type) spec builders --------

   IDS::specificationType BuildConcreteMaterialSpec(const CIdsExportOptions& options, const std::string& name,
                                                    const ConcreteMaterial& m)
   {
      IDS::applicabilityType applicability = MakeApplicability("IFCMATERIAL", nullptr, "1");
      AddAttr(applicability, "Name", SimpleValue(name));

      IDS::requirements requirements;
      const double fcR = RoundToTenthKsi_SysUnits(m.fc28);
      requirements.property().push_back(NumProperty(options, "Pset_MaterialConcrete", "CompressiveStrength",
         "IFCPRESSUREMEASURE", fcR, StressInstruction("f'c", fcR)));
      requirements.property().push_back(NumProperty(options, "Pset_MaterialConcrete", "MaxAggregateSize",
         "IFCPOSITIVELENGTHMEASURE", m.max_agg_size, LengthInstruction("max aggregate", m.max_agg_size)));

      return MakeSpec(std::move(applicability), "Concrete material - " + name, {},
         "Concrete compressive strength from PGSuper", std::move(requirements));
   }

   IDS::specificationType BuildStrandMaterialSpec(const CIdsExportOptions& options, const std::string& name,
                                                  const StrandMaterial& m)
   {
      IDS::applicabilityType applicability = MakeApplicability("IFCMATERIAL", nullptr, "1");
      AddAttr(applicability, "Name", SimpleValue(name));

      IDS::requirements requirements;
      requirements.property().push_back(NumProperty(options, "Pset_MaterialSteel", "YieldStress",
         "IFCPRESSUREMEASURE", m.fy, StressInstruction("fpy", m.fy)));
      requirements.property().push_back(NumProperty(options, "Pset_MaterialSteel", "UltimateStress",
         "IFCPRESSUREMEASURE", m.fpu, StressInstruction("fpu", m.fpu)));
      requirements.property().push_back(NumProperty(options, "Pset_MaterialSteel", "UltimateStrain",
         "IFCPOSITIVERATIOMEASURE", 0.035)); // ASTM A416, hard-coded by the IFC exporter
      requirements.property().push_back(StrProperty("Pset_MaterialSteel", "StructuralGrade", m.grade_label));
      requirements.property().push_back(NumProperty(options, "usBrPset_ACI_TendonMaterial", "TendonGrade",
         "IFCPRESSUREMEASURE", m.fpu, StressInstruction("fpu", m.fpu)));
      requirements.property().push_back(StrProperty("usBrPset_ACI_TendonMaterial", "Specification", "ASTM A416 (AASHTO M203)"));

      return MakeSpec(std::move(applicability), "Strand material - " + name, {},
         "Prestressing strand material properties from PGSuper", std::move(requirements));
   }

   IDS::specificationType BuildRebarMaterialSpec(const CIdsExportOptions& options, const std::string& name,
                                                 const RebarMaterial& m)
   {
      IDS::applicabilityType applicability = MakeApplicability("IFCMATERIAL", nullptr, "0");
      AddAttr(applicability, "Name", SimpleValue(name));

      IDS::requirements requirements;
      requirements.property().push_back(NumProperty(options, "Pset_MaterialSteel", "YieldStress",
         "IFCPRESSUREMEASURE", m.fy, StressInstruction("fy", m.fy)));
      requirements.property().push_back(NumProperty(options, "Pset_MaterialSteel", "UltimateStress",
         "IFCPRESSUREMEASURE", m.fpu, StressInstruction("fu", m.fpu)));
      requirements.property().push_back(NumProperty(options, "Pset_MaterialSteel", "UltimateStrain", "IFCPOSITIVERATIOMEASURE", m.eu));
      requirements.property().push_back(StrProperty("Pset_MaterialSteel", "StructuralGrade", name));
      requirements.property().push_back(NumProperty(options, "usBrPset_ACI_ReinforcingMaterial", "ReinforcingGrade",
         "IFCPRESSUREMEASURE", m.fy, StressInstruction("fy", m.fy)));
      requirements.property().push_back(StrProperty("usBrPset_ACI_ReinforcingMaterial", "Specification", m.spec));
      requirements.property().push_back(StrProperty("usBrPset_ACI_ReinforcingMaterial", "SpecificationVersion", m.edition));

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

   IDS::specificationType BuildBeamTypeSpec(const std::string& name)
   {
      IDS::applicabilityType applicability = MakeApplicability("IFCBEAMTYPE", nullptr, "1");
      AddAttr(applicability, "Name", SimpleValue(name));

      IDS::requirements requirements;
      requirements.property().push_back(StrProperty("Pset_ConcreteElementGeneral", "AssemblyPlace", "FACTORY"));
      requirements.property().push_back(StrProperty("Pset_ConcreteElementGeneral", "CastingMethod", "PRECAST"));

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
                  BuildBeamSpec(options, effectiveBeamId, sdv, predefinedType));

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
                        RebarSpecification(pRebar), RebarSpecificationEdition(pRebar) };
                  }
               }
            }
         }
      }

      if (options.include_beam_properties)
      {
         for (const auto& name : beamTypeNames)
            specifications.specification().push_back(BuildBeamTypeSpec(name));
      }

      if (options.include_concrete_material)
      {
         for (const auto& [name, m] : concreteMaterials)
            specifications.specification().push_back(BuildConcreteMaterialSpec(options, name, m));
      }

      if (options.include_strands)
      {
         for (const auto& [name, m] : strandMaterials)
         {
            specifications.specification().push_back(BuildStrandMaterialSpec(options, name, m));
            specifications.specification().push_back(BuildTendonTypeSpec(options, name, m));
         }
         if (anyStrandsAnywhere)
         {
            specifications.specification().push_back(BuildGlobalTendonSpec());
            specifications.specification().push_back(BuildGlobalTendonBundleSpec());
         }
      }

      if (options.include_rebar)
      {
         for (const auto& [name, m] : rebarMaterials)
            specifications.specification().push_back(BuildRebarMaterialSpec(options, name, m));
         if (anyRebarAnywhere)
         {
            specifications.specification().push_back(BuildGlobalRebarSpec());
            specifications.specification().push_back(BuildGlobalCageSpec());
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
   catch (const xml_schema::exception&)
   {
      return false;
   }
   catch (...)
   {
      return false;
   }
}
