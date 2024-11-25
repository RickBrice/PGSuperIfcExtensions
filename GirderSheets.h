///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright © 1999-2024  Washington State Department of Transportation
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

template <typename Schema>
std::vector<typename Schema::IfcDocumentReference*> GetGirderSheets()
{
   auto sheet1 = new Schema::IfcDocumentReference(
      std::string("https://wsdot.wa.gov/publications/fulltext/Bridge/Web_BSD/5.6_A4_1.PDF"), 
      boost::none, 
      std::string("WF Girder Details 1 of 5"),
      std::string("Standard Prestressed Concrete Girders"),
      nullptr);

   auto sheet2 = new Schema::IfcDocumentReference(
      std::string("https://wsdot.wa.gov/publications/fulltext/Bridge/Web_BSD/5.6_A4_2.PDF"),
      boost::none,
      std::string("WF Girder Details 2 of 5"),
      std::string("Standard Prestressed Concrete Girders"),
      nullptr);

   auto sheet3 = new Schema::IfcDocumentReference(
      std::string("https://wsdot.wa.gov/publications/fulltext/Bridge/Web_BSD/5.6_A4_3.PDF"),
      boost::none,
      std::string("WF Girder Details 3 of 5"),
      std::string("Standard Prestressed Concrete Girders"),
      nullptr);

   auto sheet4 = new Schema::IfcDocumentReference(
      std::string("https://wsdot.wa.gov/publications/fulltext/Bridge/Web_BSD/5.6_A4_4.PDF"),
      boost::none,
      std::string("WF Girder Details 4 of 5"),
      std::string("Standard Prestressed Concrete Girders"),
      nullptr);

   auto sheet5 = new Schema::IfcDocumentReference(
      std::string("https://wsdot.wa.gov/publications/fulltext/Bridge/Web_BSD/5.6_A4_5.PDF"),
      boost::none,
      std::string("WF Girder Details 5 of 5"),
      std::string("Standard Prestressed Concrete Girders"),
      nullptr);

   std::vector<typename Schema::IfcDocumentReference*> sheets({ sheet1,sheet2,sheet3,sheet4,sheet5 });
   return sheets;
}


template <typename Schema>
void AssociateDocuments(IfcHierarchyHelper<Schema>& file, typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr objects)
{
   typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr related_objects(new aggregate_of<typename Schema::IfcDefinitionSelect>());
   for (auto& object : *objects)
   {
      related_objects->push(object);
   }

   auto sheets = GetGirderSheets<Schema>();

   for (auto& sheet : sheets)
   {
      auto rel_associates_document = new Schema::IfcRelAssociatesDocument(IfcParse::IfcGlobalId(), nullptr, std::string("Standard Girder Plans"), boost::none, related_objects, sheet);
      file.addEntity(rel_associates_document);
   }
}
