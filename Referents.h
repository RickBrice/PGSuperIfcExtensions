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

#include "Properties.h"
#include <IFace\Bridge.h>
#include <PsgLib\GirderLabel.h>

#include <WBFLCogo\CogoHelpers.h>
#include <Units\StationFormat.h>

template <typename Schema>
typename Schema::IfcRelNests GetReferentNest(hierarchy_helper<Schema>& file, typename Schema::IfcAlignment alignment)
{
   auto nests = alignment.IsNestedBy();
   for (auto& nest : nests)
   {
      auto related_objects = nest.RelatedObjects();
      for (auto& related_object : related_objects)
      {
         if (auto referent = related_object.template as<IfcSchema::IfcReferent>())
         {
            return nest;
         }
      }
   }

   std::vector<typename Schema::IfcObjectDefinition> referents;
   auto rel_nests = file.create<typename Schema::IfcRelNests>().initialize(ifcopenshell::global_id(), {}, std::nullopt, std::string("Nests referents"), alignment, referents);

   return rel_nests;
}

// Creates a positioning referent for the given alignment. IfcAlignment <-> IfcRelPositions <-> IfcReferent
// The alignment positions the referent (otherwise we don't know which alignment the stationing applies to)
template <typename Schema>
typename Schema::IfcReferent CreatePositioningReferent(hierarchy_helper<Schema>& file, typename Schema::IfcAlignment alignment, std::string name, double station, typename Schema::IfcLinearPlacement placement)
{
   auto referent = file.create<typename Schema::IfcReferent>().initialize(ifcopenshell::global_id(), {}, name, std::nullopt, std::nullopt, placement, {}, Schema::IfcReferentTypeEnum::IfcReferentType_POSITION);


   // create and assign Pset_Stationing
   std::vector<typename Schema::IfcProperty> pset_station_properties;
   pset_station_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("Station"), std::nullopt, file.create<typename Schema::IfcLengthMeasure>().initialize(station), {}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("Pset_Stationing"), std::nullopt, pset_station_properties);


   std::vector<typename Schema::IfcObjectDefinition> referents;
   referents.push_back(referent);

   auto rel_defines_by_properties = file.create<typename Schema::IfcRelDefinesByProperties>().initialize(ifcopenshell::global_id(), {}, std::string("Relates pier station properties to referent"), std::nullopt, referents, property_set);


   // IfcAlignment <-> IfcRelPositions <-> IfcReferent
   // the alignment positions the referent (otherwise we don't know which alignment the stationing applies to)
   file.addRelatedObject<typename Schema::IfcRelPositions>(alignment, referent);

   return referent;
}

//
//template <typename Schema>
//void AddPositioningReferent(hierarchy_helper<Schema>& file, typename Schema::IfcAlignment alignment, typename Schema::IfcReferent referent,typename Schema::IfcProduct product)
//{
//
//   typename Schema::IfcRelNests nest = GetReferentNest<Schema>(file, alignment);
//   auto related_objects = nest->RelatedObjects();
//   related_objects.push_back(referent);
//   nest->setRelatedObjects(related_objects);
//   //std::sort(related_objects->begin(), related_objects->end(),
//   //   [](typename Schema::IfcObjectDefinition obj1, typename Schema::IfcObjectDefinition obj2)
//   //   {
//   //      typename Schema::IfcReferent ref1 = obj1->as<typename Schema::IfcReferent>();
//   //      typename Schema::IfcReferent ref2 = obj2->as<typename Schema::IfcReferent>();
//   //      if (ref1 && ref2)
//   //      {
//   //         typename Schema::IfcReal value1 = GetProperty<Schema, Schema::IfcReal>(ref1, "Pset_Stationing", "Station");
//   //         typename Schema::IfcReal value2 = GetProperty<Schema, Schema::IfcReal>(ref2, "Pset_Stationing", "Station");
//   //         if (value1 && value2)
//   //         {
//   //            return (double)(*value1) < (double)(*value2);
//   //         }
//   //      }
//   //      return false;
//   //   });
//   nest->setRelatedObjects(related_objects);
//}

// Assigns Pset_LinearReferencingMethod to the alignment, declaring that stationing (LRMName "station-point")
// is measured in feet as an absolute distance along the alignment (LRMType PEnum_LRMType.LRM_ABSOLUTE).
template <typename Schema>
void DefineLinearReferencingMethod(hierarchy_helper<Schema>& file)
{
   auto alignment = file.getSingle<typename Schema::IfcAlignment>();

   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("LRMName"), std::nullopt, file.create<typename Schema::IfcLabel>().initialize(std::string("station-point")), {}));

   // PEnum_LRMType
   std::vector<std::string> enum_values{ "LRM_ABSOLUTE","LRM_INTERPOLATIVE","LRM_RELATIVE","LRM_USERDEFINED" };
   auto property_enum_values = createPropertyEnumeration<Schema>(file, "PEnum_LRMType", enum_values); // creates an IfcPropertyEnumeration
   auto lrm_type = createPropertyEnumeratedValue<Schema>(file, "LRMType", property_enum_values, "LRM_ABSOLUTE"); // creates an IfcPropertyEnumeratedValue
   list_of_properties.push_back(lrm_type);

   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("LRMUnit"), std::nullopt, file.create<typename Schema::IfcLabel>().initialize(std::string("foot")), {}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("Pset_LinearReferencingMethod"), std::nullopt, list_of_properties);


   AddPropertySet(file, alignment, property_set);
}

template <typename Schema>
void CreateAlignmentStartStationReferent(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options)
{
   USES_CONVERSION;

   auto directrix = GetAlignmentDirectrix(file, options);

   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   auto station_format = pDisplayUnits->GetStationFormat();

   // get stationing information
   Float64 startStation, startElevation, startGrade;
   auto startPoint = GetAlignmentStartPoint(pBroker, &startStation, &startElevation, &startGrade);

   // Referents must be in order so start with the start of alignment referent

   //
   // Referent at start of alignment
   //

   // Referent position
   auto point_on_alignment = file.create<typename Schema::IfcPointByDistanceExpression>().initialize(
      file.create<typename Schema::IfcLengthMeasure>().initialize(0.0),
      std::nullopt, std::nullopt, std::nullopt,
      directrix);
   auto relative_placement = file.create<typename Schema::IfcAxis2PlacementLinear>().initialize(point_on_alignment, {}, {});
   auto referent_placement = file.create<typename Schema::IfcLinearPlacement>().initialize({}, relative_placement, {});

   // Create referent
   auto start_station_referent = file.create<typename Schema::IfcReferent>().initialize(ifcopenshell::global_id(), {}, std::string("Start of alignment station"), std::nullopt, std::nullopt, referent_placement, {}, Schema::IfcReferentTypeEnum::IfcReferentType_STATION);

   // Define properties for Pset_Stationing
   std::vector<typename Schema::IfcProperty> pset_station_properties;
   pset_station_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("Station"), std::nullopt, file.create<typename Schema::IfcLengthMeasure>().initialize(startStation), {}));

   // Create Pset and assign properties
   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("Pset_Stationing"), std::nullopt, pset_station_properties);


   // Assign the property set to the referent
   std::vector<typename Schema::IfcObjectDefinition> referents;
   referents.push_back(start_station_referent);

   auto rel_defines_by_properties = file.create<typename Schema::IfcRelDefinesByProperties>().initialize(ifcopenshell::global_id(), {}, std::string("Relates start station properties to referent"), std::nullopt, referents, property_set);


   //
   // Nest the referent to alignment
   //
   auto alignment = file.getSingle<typename Schema::IfcAlignment>();
   typename Schema::IfcRelNests nest = GetReferentNest<Schema>(file, alignment);
   auto related_objects = nest.RelatedObjects();
   related_objects.push_back(start_station_referent);
   nest.setRelatedObjects(related_objects);
}

// Returns the key point label for the referent at the start of 'segment', given the segment that
// precedes it. prev_segment is nullptr for the first segment (start of alignment) and segment is
// nullptr for the point following the last segment (end of alignment).
//
// P.C./P.T./P.I. mark a tangent-to-curve/curve-to-tangent/tangent-to-tangent point with no spiral.
// T.S./S.T./S.C./C.S. mark the tangent/spiral and spiral/curve transitions of a spiraled curve.
// P.C.C. marks a compound curve (two circular arcs meeting directly, no intervening tangent).
template <typename Schema>
std::string GetHorizontalKeyPointLabel(typename Schema::IfcAlignmentSegment prev_segment, typename Schema::IfcAlignmentSegment segment)
{
   if (!prev_segment) return "P.O.B."; // Point of Beginning
   if (!segment) return "P.O.E.";      // Point of Ending

   auto prev_type = prev_segment.DesignParameters().template as<typename Schema::IfcAlignmentHorizontalSegment>().PredefinedType();
   auto type = segment.DesignParameters().template as<typename Schema::IfcAlignmentHorizontalSegment>().PredefinedType();

   if (prev_type == Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_LINE)
   {
      if (type == Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_LINE)        return "P.I.";
      if (type == Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CIRCULARARC) return "P.C.";
      if (type == Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID)    return "T.S.";
   }
   else if (prev_type == Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CIRCULARARC)
   {
      if (type == Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_LINE)        return "P.T.";
      if (type == Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CIRCULARARC) return "P.C.C.";
      if (type == Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID)    return "C.S.";
   }
   else if (prev_type == Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID)
   {
      if (type == Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_LINE)        return "S.T.";
      if (type == Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CIRCULARARC) return "S.C.";
   }

   CHECK(false); // unexpected horizontal segment transition
   return "";
}

// Returns the key point label for the referent at the start of 'segment', given the segment that
// precedes it. Only called for interior transitions (both prev_segment and segment are non-null) --
// the vertical profile's start/end points coincide with P.O.B./P.O.E. of the horizontal alignment
// and are not labeled again here.
//
// B.V.C./E.V.C. mark the beginning/end of a parabolic vertical curve. P.V.I. marks a grade break
// with no vertical curve. V.C.C. marks a compound vertical curve (two parabolic arcs meeting
// directly at their PVI, no intervening constant grade).
template <typename Schema>
std::string GetVerticalKeyPointLabel(typename Schema::IfcAlignmentSegment prev_segment, typename Schema::IfcAlignmentSegment segment)
{
   auto prev_type = prev_segment.DesignParameters().template as<typename Schema::IfcAlignmentVerticalSegment>().PredefinedType();
   auto type = segment.DesignParameters().template as<typename Schema::IfcAlignmentVerticalSegment>().PredefinedType();

   if (prev_type == Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CONSTANTGRADIENT)
   {
      if (type == Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CONSTANTGRADIENT) return "P.V.I.";
      if (type == Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_PARABOLICARC)     return "B.V.C.";
   }
   else if (prev_type == Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_PARABOLICARC)
   {
      if (type == Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CONSTANTGRADIENT) return "E.V.C.";
      if (type == Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_PARABOLICARC)     return "V.C.C.";
   }

   CHECK(false); // unexpected vertical segment transition
   return "";
}

// Creates an IfcReferent at 'distance_along' the alignment directrix, named "<label> (<station>)", with
// its Pset_Stationing.Station property set. Does not nest the referent to anything -- the caller does that.
template <typename Schema>
typename Schema::IfcReferent CreateKeyPointReferent(hierarchy_helper<Schema>& file, typename Schema::IfcCurve directrix, const std::string& label, double distance_along, double station, const WBFL::Units::StationFormat& station_format)
{
   USES_CONVERSION;

   auto point_on_alignment = file.create<typename Schema::IfcPointByDistanceExpression>().initialize(
      file.create<typename Schema::IfcLengthMeasure>().initialize(distance_along),
      std::nullopt, std::nullopt, std::nullopt,
      directrix);
   auto relative_placement = file.create<typename Schema::IfcAxis2PlacementLinear>().initialize(point_on_alignment, {}, {});
   auto referent_placement = file.create<typename Schema::IfcLinearPlacement>().initialize({}, relative_placement, {});

   std::ostringstream os;
   os << label << " (" << T2A(WBFL::COGO::Station(station).AsString(station_format).c_str()) << ")";

   auto referent = file.create<typename Schema::IfcReferent>().initialize(ifcopenshell::global_id(), {}, os.str(), std::nullopt, std::nullopt, referent_placement, {}, Schema::IfcReferentTypeEnum::IfcReferentType_POSITION);


   std::vector<typename Schema::IfcProperty> pset_station_properties;
   pset_station_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("Station"), std::nullopt, file.create<typename Schema::IfcLengthMeasure>().initialize(station), {}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("Pset_Stationing"), std::nullopt, pset_station_properties);


   std::vector<typename Schema::IfcObjectDefinition> referents;
   referents.push_back(referent);

   auto rel_defines_by_properties = file.create<typename Schema::IfcRelDefinesByProperties>().initialize(ifcopenshell::global_id(), {}, std::string("Relates key point station properties to referent"), std::nullopt, referents, property_set);


   return referent;
}

// Creates IfcReferent key-point markers for every segment transition in the horizontal alignment and
// vertical profile (P.O.B., P.O.E., P.C., P.T., P.I., T.S., S.T., S.C., C.S., P.C.C. for horizontal;
// B.V.C., E.V.C., P.V.I., V.C.C. for vertical) and nests them with two new IfcRelNests: the horizontal
// key points to 'horizontal', and the vertical key points to 'vertical'.
//
// horizontal/nests_horizontal_segments and vertical/nests_vertical_segments are the layouts and
// segment nests created by CreateHorizontalAlignment/CreateVerticalProfile. The last entry in each
// segment nest is the zero-length terminator segment (CT 4.1.7.1.1.2) and is skipped.
template <typename Schema>
void UpdateKeyPointReferents(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options,
   typename Schema::IfcAlignmentHorizontal horizontal, typename Schema::IfcRelNests nests_horizontal_segments,
   typename Schema::IfcAlignmentVertical vertical, typename Schema::IfcRelNests nests_vertical_segments)
{
   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   auto station_format = pDisplayUnits->GetStationFormat();

   Float64 startStation, startElevation, startGrade;
   GetAlignmentStartPoint(pBroker, &startStation, &startElevation, &startGrade);

   auto directrix = GetAlignmentDirectrix<Schema>(file, options);

   std::vector<typename Schema::IfcObjectDefinition> new_horizontal_referents;
   std::vector<typename Schema::IfcObjectDefinition> new_vertical_referents;

   // horizontal key points
   {
      auto segments = nests_horizontal_segments.RelatedObjects();
      auto end = std::prev(segments.end()); // skip the zero-length terminator segment

      double distance_along = 0.0;
      typename Schema::IfcAlignmentSegment prev_segment;
      for (auto iter = segments.begin(); iter != end; iter++)
      {
         auto segment = iter->template as<typename Schema::IfcAlignmentSegment>();
         auto dp = segment.DesignParameters().template as<typename Schema::IfcAlignmentHorizontalSegment>();

         auto label = GetHorizontalKeyPointLabel<Schema>(prev_segment, segment);
         new_horizontal_referents.push_back(CreateKeyPointReferent<Schema>(file, directrix, label, distance_along, startStation + distance_along, station_format));

         distance_along += dp.SegmentLength();
         prev_segment = segment;
      }

      auto label = GetHorizontalKeyPointLabel<Schema>(prev_segment, {});
      new_horizontal_referents.push_back(CreateKeyPointReferent<Schema>(file, directrix, label, distance_along, startStation + distance_along, station_format));
   }

   // vertical key points (interior transitions only -- see GetVerticalKeyPointLabel)
   {
      auto segments = nests_vertical_segments.RelatedObjects();
      auto end = std::prev(segments.end()); // skip the zero-length terminator segment

      typename Schema::IfcAlignmentSegment prev_segment;
      for (auto iter = segments.begin(); iter != end; iter++)
      {
         auto segment = iter->template as<typename Schema::IfcAlignmentSegment>();
         if (prev_segment)
         {
            auto dp = segment.DesignParameters().template as<typename Schema::IfcAlignmentVerticalSegment>();
            auto label = GetVerticalKeyPointLabel<Schema>(prev_segment, segment);
            new_vertical_referents.push_back(CreateKeyPointReferent<Schema>(file, directrix, label, dp.StartDistAlong(), startStation + dp.StartDistAlong(), station_format));
         }
         prev_segment = segment;
      }
   }

   auto alignment = file.getSingle<typename Schema::IfcAlignment>();

   auto nests_horizontal_referents = file.create<typename Schema::IfcRelNests>().initialize(ifcopenshell::global_id(), {}, std::nullopt, std::string("Nests horizontal key point referents with alignment"), alignment, new_horizontal_referents);


   auto nests_vertical_referents = file.create<typename Schema::IfcRelNests>().initialize(ifcopenshell::global_id(), {}, std::nullopt, std::string("Nests vertical key point referents with alignment"), alignment, new_vertical_referents);

}
