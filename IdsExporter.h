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

// Reusable, UI-independent generator for a project-specific buildingSMART IDS
// (Information Delivery Specification) describing the girder concrete strength design
// values PGSuper computed. Takes only a broker, an options struct, and an output sink,
// so it can be driven from the standalone IDS Data Exporter plug-in today and from the
// IFC exporter (co-emitting a .ids next to the .ifc) in a future revision.
//
// This header deliberately exposes no xsd-cxx / Xerces types; the generated IDS schema
// binding is an implementation detail of IdsExporter.cpp.

#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <PsgLib\Keys.h>

namespace WBFL { namespace EAF { class Broker; } }

class CIdsExportOptions
{
public:
   // Master on/off switch for the IDS tab of the export dialog. IDS export only ever
   // runs alongside an IFC build (see CExportOptionsSheet / PGSuperDataExporter.cpp),
   // so this is the single gate for whether that second file gets written at all.
   bool enabled = false;

   // How a pinned numeric requirement constrains the value found in the IFC model.
   // Applies to every exact-pinned value the generator emits (concrete strengths,
   // jacking force/stress, material strengths, geometry, cambers), not just strengths.
   enum class StrengthMatch
   {
      Exact,   // <ids:simpleValue> with the design value (validator applies its tolerance)
      Minimum  // <xs:restriction base="xs:double"><xs:minInclusive .../></xs:restriction>
   };

   // How each <specification> pins to one specific IfcBeam.
   enum class BeamId
   {
      Name,    // <ids:attribute> Name == "Span i, Girder X"  (the IfcBeam.Name string)
      GlobalId // <ids:attribute> GlobalId == the IfcBeam GlobalId
   };

   StrengthMatch strength_match = StrengthMatch::Exact;
   BeamId        beam_id        = BeamId::Name;

   bool include_release_strength        = true; // Pset_PrecastConcreteElementGeneral.ReleaseStrength
   bool include_transportation_strength = true; // Pset_PrecastConcreteElementGeneral.TransportationStrength

   // Categories of the prestressed-girder export the IDS should assert. Each maps to
   // a group of <specification>s. See IdsExporter.cpp for exactly what each emits.
   bool include_beam_properties  = true; // Pset_BeamCommon, Pset_(Precast)ConcreteElementGeneral,
                                         // usBrPset_PrecastConcreteBeam, IfcBeamType, classification, material
   bool include_quantities       = true; // Qto_BeamBaseQuantities (presence)
   bool include_concrete_material = true; // concrete IfcMaterial + Pset_MaterialConcrete
   bool include_strands          = true; // IfcTendon / IfcTendonType / strand IfcMaterial / Strands assembly
   bool include_rebar            = true; // IfcReinforcingBar / rebar IfcMaterial / reinforcement cage

   // Consulted only when beam_id == BeamId::GlobalId. Filled in by PGSuperDataExporter.cpp
   // from the IfcBeam GlobalIds the same CIfcExporter::BuildModel call just created.
   std::map<CSegmentKey, std::string> global_id_by_segment;

   // <ids:info> metadata. Empty fields are given sensible defaults by the worker
   // (title <- bridge name; date <- today). "author" is emitted only if it looks like
   // an e-mail address (the IDS schema rejects anything else).
   CString title;
   CString copyright;
   CString version;
   CString description;
   CString author;
   CString purpose;
   CString milestone;

   // Prefix for each generated girder <specification> name, e.g. "Girder design values".
   CString specification_name_prefix = _T("Girder design values");
};

class CIdsExporter
{
public:
   CIdsExporter() = default;
   ~CIdsExporter() = default;

   // Opens strFilePath and delegates to the ostream overload. Returns false on any error.
   bool BuildSpecification(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIdsExportOptions& options, const CString& strFilePath);

   // Writes the IDS document to os. Returns false on any error (including analysis
   // errors raised by the broker while gathering strengths).
   bool BuildSpecification(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIdsExportOptions& options, std::ostream& os);
};
