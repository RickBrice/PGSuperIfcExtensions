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

// Helpers for building buildingSMART IDS documents with the xsd-cxx binding (Schema/ids-binding.hxx),
// shared by the design-value IDS (IdsExporter.cpp) and the general IDS written from a mapping table
// (IfcTableIds.cpp). Include only in .cpp files: the binding is an implementation detail.

#include "Schema/ids-binding.hxx"

#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>

#include <atlconv.h>
#include <ctime>
#include <optional>
#include <string>
#include <vector>

namespace ids_builder
{
   inline const char* const IDS_NAMESPACE = "http://standards.buildingsmart.org/IDS";
   inline const char* const IDS_SCHEMA_LOCATION = "http://standards.buildingsmart.org/IDS/1.0/ids.xsd";
   inline const char* const XS_NAMESPACE = "http://www.w3.org/2001/XMLSchema";

   inline std::string ToUtf8(LPCTSTR s)
   {
      if (s == nullptr) return std::string();
      CW2A conv(s, CP_UTF8);
      return std::string((LPCSTR)conv);
   }

   inline std::string ToUtf8(const CString& s) { return ToUtf8((LPCTSTR)s); }

   inline xml_schema::date Today()
   {
      std::time_t t = std::time(nullptr);
      std::tm lt{};
      localtime_s(&lt, &t);
      return xml_schema::date(lt.tm_year + 1900, static_cast<unsigned short>(lt.tm_mon + 1), static_cast<unsigned short>(lt.tm_mday));
   }

   inline IDS::idsValue SimpleValue(const std::string& text)
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

   // <ids:value><xs:restriction base="..."><xs:<facet> value="..."/>...</xs:restriction></ids:value>, one facet element per value
   inline IDS::idsValue RestrictionValue(const char* base, const char* facet, const std::vector<std::string>& facetValues)
   {
      IDS::idsValue v;
      xercesc::DOMDocument& doc = v.dom_document();

      XStr xsNs(XS_NAMESPACE);

      xercesc::DOMElement* restriction = doc.createElementNS(xsNs, XStr("xs:restriction"));
      restriction->setAttribute(XStr("base"), XStr(base));

      for (const auto& facetValue : facetValues)
      {
         xercesc::DOMElement* facetElem = doc.createElementNS(xsNs, XStr(facet));
         facetElem->setAttribute(XStr("value"), XStr(facetValue.c_str()));
         restriction->appendChild(facetElem);
      }

      v.any(restriction);
      return v;
   }

   inline IDS::idsValue RestrictionValue(const char* base, const char* facet, const std::string& facetValue)
   {
      return RestrictionValue(base, facet, std::vector<std::string>{ facetValue });
   }

   // one of the values: <xs:restriction base="xs:string"><xs:enumeration value="..."/>...
   inline IDS::idsValue EnumerationValue(const std::vector<std::string>& values)
   {
      return RestrictionValue("xs:string", "xs:enumeration", values);
   }

   enum class Card { Required, Optional };

   inline const char* CardText(Card c) { return c == Card::Optional ? "optional" : "required"; }

   inline IDS::property PropertyReq(const std::string& pset, const std::string& baseName, const char* ifcDataType,
                                    std::optional<IDS::idsValue> value, Card card, const std::string& instruction = {}, const std::string& uri = {})
   {
      IDS::property p(SimpleValue(pset), SimpleValue(baseName));
      if (ifcDataType) p.dataType(IDS::upperCaseName(ifcDataType));
      p.cardinality(IDS::conditionalCardinality(CardText(card)));
      if (!instruction.empty()) p.instructions(instruction);
      if (!uri.empty()) p.uri(uri);
      if (value) p.value(std::move(*value));
      return p;
   }

   inline IDS::classification ClassificationReq(const std::string& system, const std::string& code, const std::string& uri = {})
   {
      IDS::classification c(SimpleValue(system));
      c.value(SimpleValue(code));           // matches IfcClassificationReference.Identification
      c.cardinality(IDS::conditionalCardinality("required"));
      if (!uri.empty()) c.uri(uri);
      return c;
   }

   inline IDS::material MaterialReq(const std::string& nameOrCategory)
   {
      IDS::material m;
      m.value(SimpleValue(nameOrCategory)); // matched against IfcMaterial.Name and .Category
      m.cardinality(IDS::conditionalCardinality("required"));
      return m;
   }

   inline IDS::attribute AttributeReq(const std::string& name, std::optional<IDS::idsValue> value, const std::string& instruction = {})
   {
      IDS::attribute a(SimpleValue(name));
      a.cardinality(IDS::conditionalCardinality("required"));
      if (!instruction.empty()) a.instructions(instruction);
      if (value) a.value(std::move(*value));
      return a;
   }

   inline IDS::applicabilityType MakeApplicability(const char* ifcClass, const char* predefinedType, const std::string& minOccurs)
   {
      IDS::entityType entity(SimpleValue(ifcClass));
      if (predefinedType) entity.predefinedType(SimpleValue(predefinedType));

      IDS::applicabilityType applicability;
      applicability.entity(entity);
      applicability.minOccurs(minOccurs);
      applicability.maxOccurs(std::string("unbounded"));
      return applicability;
   }

   inline void AddAttr(IDS::applicabilityType& applicability, const std::string& name, IDS::idsValue value)
   {
      IDS::attributeType attribute(SimpleValue(name));
      attribute.value(std::move(value));
      applicability.attribute().push_back(attribute);
   }

   inline IDS::specificationType MakeSpec(IDS::applicabilityType applicability, const std::string& name,
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

   // writes the document with the ids and xs namespace prefixes
   inline void Write(std::ostream& os, const IDS::ids& document)
   {
      xml_schema::namespace_infomap map;
      map["ids"].name = IDS_NAMESPACE;
      map["ids"].schema = IDS_SCHEMA_LOCATION;
      map["xs"].name = XS_NAMESPACE;
      IDS::ids_(os, document, map, "UTF-8");
   }
}
