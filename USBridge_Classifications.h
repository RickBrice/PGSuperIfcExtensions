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

#include <IFace/Tools.h>
#include <IFace/PointOfInterest.h>
#include <IFace\Bridge.h>
#include <IFace\AnalysisResults.h>
#include <IFace\Intervals.h>
#include "Units.h"
#include "IfcExporter.h"
#include "Properties.h"
#include "bSDD.h"

/// @brief Checks if an object has the specified classification
/// @tparam Schema 
/// @param object 
/// @param identifier 
/// @return 
template <typename Schema>
bool HasClassification(typename Schema::IfcObjectDefinition object, std::string identifier)
{
   auto associations = object.HasAssociations();
   for (auto& rel : associations)
   {
      auto rel_associates_classification = rel.template as<typename Schema::IfcRelAssociatesClassification>();
      if (rel_associates_classification)
      {
         auto classification_reference = rel_associates_classification.RelatingClassification().template as<typename Schema::IfcClassificationReference>();
         if (classification_reference && classification_reference.Identification().value_or("") == identifier)
            return true;
      }
   }

   return false;
}


template <typename Schema>
void AssociateClassification(hierarchy_helper<Schema>& file, typename Schema::IfcClassificationSelect classification, typename Schema::IfcDefinitionSelect related_object)
{
   auto associations = related_object.template as<typename Schema::IfcObjectDefinition>().HasAssociations();
   for (auto& association : associations)
   {
      auto rel_classification = association.template as<typename Schema::IfcRelAssociatesClassification>();
      if (rel_classification && rel_classification.RelatingClassification() == classification)
      {
         auto related_objects = rel_classification.RelatedObjects();
         related_objects.push_back(related_object);
         rel_classification.setRelatedObjects(related_objects);
         return;
      }
   }

   // if we get this far, there was not already a classicification association for this product, so we will create a new one
   std::vector<typename Schema::IfcDefinitionSelect> related_objects;
   related_objects.push_back(related_object);
   auto rel_classification = file.create<typename Schema::IfcRelAssociatesClassification>();
   rel_classification.setGlobalId(ifcopenshell::global_id());
   rel_classification.setRelatedObjects(related_objects);
   rel_classification.setRelatingClassification(classification);
}

template <typename Schema>
void Classify_usBridge_ReinforcingBarType(hierarchy_helper<Schema>& file, typename Schema::IfcReinforcingBarType rebar_type)
{
   // usBridge does not define a classification for IfcReinforcingBarType.
   // Classification it is unclear if classification association is inherited from the type to the instance.
   // When usBridge defines one, declare it in the mapping table and use Classify (IfcPropertyWriter.h)
}