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

#include "Materials.h"
#include "BarShape.h"
#include "PropertySets.h"
#include "QuantitySets.h"
#include "RebarRelationshipBatch.h"

#include <array>
#include <map>
#include <sstream>

inline double getMinBendRadius(const WBFL::Materials::Rebar* pRebar,bool bStirrup)
{
   auto radius = 0.0;
   auto db = pRebar->GetNominalDimension();
   if (pRebar->GetSize() <= WBFL::Materials::Rebar::Size::bs5 && bStirrup)
      radius = 4 * db;
   else if (pRebar->GetSize() <= WBFL::Materials::Rebar::Size::bs8)
      radius = 6 * db;
   else if (pRebar->GetSize() <= WBFL::Materials::Rebar::Size::bs11)
      radius = 8 * db;
   else
      radius = 10 * db;

   return radius;
}

// The directrix of a bar: its centerline path as an IfcIndexedPolyCurve
template <typename Schema>
typename Schema::IfcIndexedPolyCurve CreateBarDirectrix(hierarchy_helper<Schema>& file, const BarPath& path)
{
   std::vector<std::vector<double>> points;
   for (const auto& p : path.GetPoints())
      points.push_back({ p[0], p[1], p[2] });

   // IFC indices are 1-based
   std::vector<typename Schema::IfcSegmentIndexSelect> segments;
   for (const auto& s : path.GetSegments())
   {
      if (s.IsArc())
         segments.push_back(file.create<typename Schema::IfcArcIndex>().initialize({ (int)s.a + 1, (int)s.b + 1, (int)s.c + 1 }));
      else
         segments.push_back(file.create<typename Schema::IfcLineIndex>().initialize({ (int)s.a + 1, (int)s.b + 1 }));
   }

   return file.create<typename Schema::IfcIndexedPolyCurve>().initialize(file.create<typename Schema::IfcCartesianPointList3D>().initialize(points, std::nullopt), segments, std::nullopt);
}

// The IfcReinforcingBarType of a bar (see BarShape.h). It's created for the first bar with its definition and shared
// by the bars after it. Everything about the type comes from the definition, so a new bar shape only needs its
// BarDefinition:
// - the geometry is the centerline path swept by the bar diameter, mapped by the occurrences (see MapBarType)
// - BarLength is the centerline length. NominalDiameter and CrossSectionArea are the bar's
// - the mapping table's type property sets of the rebar role, with the bar in the export context (bar mark, size,
//   mass, use, and covers), and usBrPset_ACI_BarShape from the bend shape
// - the material and the classification of the rebar role
// usBridge has a bar's data on its type, not on the bar (usBrSpec_Classes.ids), so bars with different paths or
// type property values get different types. A second type with the same name gets " (2)", and so on.
template <typename Schema>
typename Schema::IfcReinforcingBarType GetBarType(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, const BarDefinition& bar,
   typename Schema::IfcReinforcingBarTypeEnum::Value predefined_type, typename Schema::IfcMaterial material, RebarRelationshipBatch<Schema>& rebar_batch)
{
   ExportContext context{ pBroker, &options };
   context.rebar = bar.rebar;
   context.bar = &bar;

   std::ostringstream key;
   key << bar.name << "|" << (int)bar.rebar->GetType() << "," << (int)bar.rebar->GetGrade() << "," << (int)bar.rebar->GetSize() << "|" << (int)predefined_type << "|" << material.id()
      << "|" << bar.path.GetSignature() << "|" << bar.shape_code << "," << bar.bend_radius;
   std::map<std::string, Float64> dimensions(bar.shape_dimensions.begin(), bar.shape_dimensions.end()); // sorted
   for (const auto& [name, value] : dimensions)
      key << "," << name << "=" << value;
   key << "|" << GetTypePropertySignature(ElementKind::Rebar, context);

   auto& session = CIfcExportSession::Current();
   if (auto found = session.bar_types.find(key.str()); found != session.bar_types.end())
      return file.instance_by_id(found->second).template as<typename Schema::IfcReinforcingBarType>();

   auto count = ++session.type_names[bar.name];
   std::string name = (count == 1 ? bar.name : bar.name + " (" + std::to_string(count) + ")");

   // geometry
   auto geometric_representation_context = file.getRepresentationContext(std::string("Model"));
   auto swept_disk_solid = file.create<typename Schema::IfcSweptDiskSolid>().initialize(CreateBarDirectrix<Schema>(file, bar.path), bar.rebar->GetNominalDimension() / 2, std::nullopt, std::nullopt, std::nullopt);
   std::vector<typename Schema::IfcRepresentationItem> representation_items{ swept_disk_solid };
   auto shape_representation = file.create<typename Schema::IfcShapeRepresentation>().initialize(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);

   // style the mapped representation - every IfcMappedItem of this bar type gets its color from it
   // main bars use the Rebar style, all others the Stirrup style (as the Rebar and Stirrup materials)
   auto surface_style = (predefined_type == Schema::IfcReinforcingBarTypeEnum::IfcReinforcingBarType_MAIN ? GetSurfaceStyle<Schema>(file, "Rebar", REBAR_COLOR) : GetSurfaceStyle<Schema>(file, "Stirrup", STIRRUP_COLOR));
   StyleShapeRepresentation<Schema>(file, shape_representation, surface_style);

   auto representation_map = file.create<typename Schema::IfcRepresentationMap>().initialize(file.addPlacement3d(), shape_representation);
   std::vector<typename Schema::IfcRepresentationMap> representation_maps{ representation_map };

   // property sets
   auto property_sets = CreateTypePropertySets<Schema>(file, ElementKind::Rebar, context);
   if (options.classify && !bar.shape_code.empty())
   {
      if (auto bar_shape = Create_usBrPset_ACI_BarShape(file, pBroker, options, bar.shape_code, bar.bend_radius, bar.shape_dimensions))
         property_sets.push_back(bar_shape);
   }
   std::optional<std::vector<typename Schema::IfcPropertySetDefinition>> has_property_sets;
   if (!property_sets.empty())
      has_property_sets = property_sets;

   auto rebar_type = file.create<typename Schema::IfcReinforcingBarType>().initialize(
      ifcopenshell::global_id(),
      {},
      name, /*Name*/
      std::nullopt, /*Description*/
      std::string("IfcReinforcingBar"), /*ApplicableOccurrence*/
      has_property_sets, /*HasPropertySets*/
      representation_maps, /*RepresentationMaps*/
      std::nullopt, /*Tag*/
      std::nullopt, /*ElementType*/
      predefined_type, /*PredefinedType*/
      bar.rebar->GetNominalDimension(), /*NominalDiameter*/
      bar.rebar->GetNominalArea(), /*CrossSectionArea*/
      bar.GetLength(), /*BarLength*/
      std::nullopt, /*BarSurface*/
      std::nullopt, /*BendingShapeCode*/
      std::nullopt /*BendingParameters*/
   );

   rebar_batch.Material(material, rebar_type);
   if (options.classify)
      rebar_batch.Classify(GetClassificationReferences<Schema>(file, ElementKind::Rebar, options), rebar_type);

   file.addRelatedObject<typename Schema::IfcRelDeclares>(file.getSingle<typename Schema::IfcProject>(), rebar_type);

   session.bar_types.emplace(key.str(), rebar_type.id());
   return rebar_type;
}

// The geometry of one bar: its bar type's geometry, with the type's x and y axes along x_dir and y_dir, at origin
// (in the coordinates of the occurrence's placement)
template <typename Schema>
typename Schema::IfcMappedItem MapBarType(hierarchy_helper<Schema>& file, typename Schema::IfcReinforcingBarType rebar_type, const std::array<double, 3>& origin, const std::array<double, 3>& x_dir, const std::array<double, 3>& y_dir)
{
   auto mapping_target = file.create<typename Schema::IfcCartesianTransformationOperator3D>().initialize(
      file.create<typename Schema::IfcDirection>().initialize({ x_dir[0], x_dir[1], x_dir[2] }),
      file.create<typename Schema::IfcDirection>().initialize({ y_dir[0], y_dir[1], y_dir[2] }),
      file.create<typename Schema::IfcCartesianPoint>().initialize({ origin[0], origin[1], origin[2] }),
      std::nullopt, {});
   return file.create<typename Schema::IfcMappedItem>().initialize(rebar_type.RepresentationMaps()->front(), mapping_target);
}

// A bar to write: its type, its geometry, and its name
template <typename Schema>
struct PlacedBar
{
   typename Schema::IfcReinforcingBarType type;
   typename Schema::IfcMappedItem geometry;
   std::string name;
};

// Writes the IfcReinforcingBar occurrences of a group of bars (e.g. a row or a stirrup zone): one per bar (ACI 131
// Article 7.3, a 1:1 correspondence between a real bar and an IfcReinforcingBar), or, with
// CIfcExportOptions::RebarRepresentation::Mapped, one per bar type with the geometry of all of its bars and a Count
// quantity (group_name, with the type name when the group has more than one type). An occurrence has its geometry,
// its type, and its assembly: usBridge has the bar's data on its type
template <typename Schema>
void WriteBars(hierarchy_helper<Schema>& file, const CIfcExportOptions& options, typename Schema::IfcObjectPlacement placement, typename Schema::IfcElementAssembly rebar_assembly, RebarRelationshipBatch<Schema>& rebar_batch,
   const std::vector<PlacedBar<Schema>>& bars, const std::string& group_name)
{
   auto geometric_representation_context = file.getRepresentationContext(std::string("Model"));
   auto write = [&](typename Schema::IfcReinforcingBarType type, const std::vector<typename Schema::IfcRepresentationItem>& items, const std::string& name)
      {
         auto shape_representation = file.create<typename Schema::IfcShapeRepresentation>().initialize(geometric_representation_context, std::string("Body"), std::string("MappedRepresentation"), items);
         std::vector<typename Schema::IfcRepresentation> shape_representation_list{ shape_representation };
         auto product_definition_shape = file.create<typename Schema::IfcProductDefinitionShape>().initialize(std::nullopt, std::nullopt, shape_representation_list);

         // the attributes about the bar are deprecated (they're on the type)
         auto rebar = file.create<typename Schema::IfcReinforcingBar>().initialize(ifcopenshell::global_id(), {}, name, std::nullopt, std::nullopt, placement, product_definition_shape, std::nullopt,
            std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt);

         rebar_batch.Type(type, rebar);
         rebar_batch.Aggregate(rebar_assembly, rebar);
         return rebar;
      };

   if (options.rebar_representation != CIfcExportOptions::RebarRepresentation::Mapped)
   {
      for (const auto& bar : bars)
         write(bar.type, { bar.geometry }, bar.name);
      return;
   }

   // the bars of each type, in the order of the bars
   std::vector<std::pair<typename Schema::IfcReinforcingBarType, std::vector<typename Schema::IfcRepresentationItem>>> groups;
   for (const auto& bar : bars)
   {
      auto group = std::find_if(groups.begin(), groups.end(), [&bar](const auto& g) { return g.first.id() == bar.type.id(); });
      if (group == groups.end())
         groups.emplace_back(bar.type, std::vector<typename Schema::IfcRepresentationItem>{ bar.geometry });
      else
         group->second.push_back(bar.geometry);
   }

   for (const auto& [type, items] : groups)
   {
      std::string name = (groups.size() == 1 ? group_name : group_name + ", " + std::string(type.Name().value_or("")));
      auto rebar = write(type, items, name);

      // A grouped occurrence needs its count. usBridge doesn't allow quantities on a bar, so Mapped isn't
      // usBridge conformant for this quantity (Individual is)
      if (options.classify)
         AddQto(file, rebar, Create_Qto_ReinforcingElementGroupQuantities<Schema>(file, items.size()));
   }
}
