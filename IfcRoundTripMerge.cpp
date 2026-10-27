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
#include "IfcRoundTripMerge.h"
#include "IfcImportUnits.h"
#include "IfcMappingTable.h"
#include "IfcTargetReader.h"

#include <Eigen/Dense>

#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>

namespace
{
   using UnitType = CIfcImportUnits::UnitType;

   // the unit type of a measure (IfcValue) type, or nullopt if its value has no unit
   std::optional<UnitType> unit_type_of_measure(const std::string& type)
   {
      static const std::map<std::string, UnitType> types{
         {"IfcLengthMeasure", IfcSchema::IfcUnitEnum::IfcUnit_LENGTHUNIT},
         {"IfcPositiveLengthMeasure", IfcSchema::IfcUnitEnum::IfcUnit_LENGTHUNIT},
         {"IfcNonNegativeLengthMeasure", IfcSchema::IfcUnitEnum::IfcUnit_LENGTHUNIT},
         {"IfcAreaMeasure", IfcSchema::IfcUnitEnum::IfcUnit_AREAUNIT},
         {"IfcVolumeMeasure", IfcSchema::IfcUnitEnum::IfcUnit_VOLUMEUNIT},
         {"IfcMassMeasure", IfcSchema::IfcUnitEnum::IfcUnit_MASSUNIT},
         {"IfcForceMeasure", IfcSchema::IfcUnitEnum::IfcUnit_FORCEUNIT},
         {"IfcPressureMeasure", IfcSchema::IfcUnitEnum::IfcUnit_PRESSUREUNIT},
         {"IfcPlaneAngleMeasure", IfcSchema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT},
         {"IfcPositivePlaneAngleMeasure", IfcSchema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT},
      };
      auto found = types.find(type);
      return found == types.end() ? std::nullopt : std::optional<UnitType>(found->second);
   }

   // the value attribute and unit type of a simple quantity
   std::optional<std::pair<std::string, UnitType>> quantity_value_of(const std::string& type)
   {
      if (type == "IfcQuantityLength") return std::make_pair(std::string("LengthValue"), IfcSchema::IfcUnitEnum::IfcUnit_LENGTHUNIT);
      if (type == "IfcQuantityArea") return std::make_pair(std::string("AreaValue"), IfcSchema::IfcUnitEnum::IfcUnit_AREAUNIT);
      if (type == "IfcQuantityVolume") return std::make_pair(std::string("VolumeValue"), IfcSchema::IfcUnitEnum::IfcUnit_VOLUMEUNIT);
      if (type == "IfcQuantityWeight") return std::make_pair(std::string("WeightValue"), IfcSchema::IfcUnitEnum::IfcUnit_MASSUNIT);
      return std::nullopt;
   }

   // A value as text for the log and for comparisons
   std::string value_text(const express::base& value)
   {
      if (!value)
         return "(none)";
      auto attribute = value.get_attribute_value(0);
      if (attribute.isNull())
         return "(none)";
      std::ostringstream os;
      switch (attribute.type())
      {
      case ifcopenshell::Argument_INT: os << (int64_t)attribute; break;
      case ifcopenshell::Argument_DOUBLE: os << std::setprecision(10) << (double)attribute; break;
      case ifcopenshell::Argument_BOOL: os << ((bool)attribute ? "true" : "false"); break;
      case ifcopenshell::Argument_LOGICAL:
      {
         boost::logic::tribool b = attribute;
         os << (boost::logic::indeterminate(b) ? "unknown" : (bool)b ? "true" : "false");
         break;
      }
      case ifcopenshell::Argument_STRING: os << "\"" << (std::string)attribute << "\""; break;
      default: os << "(value)"; break;
      }
      return os.str();
   }

   std::optional<double> number_of(const express::base& value)
   {
      if (!value)
         return std::nullopt;
      auto attribute = value.get_attribute_value(0);
      if (attribute.type() == ifcopenshell::Argument_DOUBLE)
         return (double)attribute;
      if (attribute.type() == ifcopenshell::Argument_INT)
         return (double)(int64_t)attribute;
      return std::nullopt;
   }

   // The registry key's element role, for the mapping table's property declarations
   std::vector<ElementKind> roles_of_key(const std::string& key)
   {
      auto prefix = key.substr(0, key.find(':'));
      if (prefix == "project") return { ElementKind::Project };
      if (prefix == "site") return { ElementKind::Site };
      if (prefix == "bridge") return { ElementKind::Bridge };
      if (prefix == "superstructure") return { ElementKind::Superstructure };
      if (prefix == "substructure") return { ElementKind::Substructure };
      if (prefix == "deck_part") return { ElementKind::DeckPart };
      if (prefix == "deck") return { ElementKind::Deck };
      if (prefix == "barrier") return { ElementKind::Barrier };
      if (prefix == "pier") return { ElementKind::Pier, ElementKind::Abutment };
      if (prefix == "girder") return { ElementKind::Girder };
      if (prefix == "bearings") return { ElementKind::Bearing };
      if (prefix == "diaphragms") return { ElementKind::Diaphragm };
      if (prefix == "tendons") return { ElementKind::TendonBundle };
      if (prefix == "bars") return { ElementKind::ReinforcementCage };
      return {};
   }

   std::string id_of_key(const std::string& key)
   {
      auto colon = key.find(':');
      return colon == std::string::npos ? std::string() : key.substr(colon + 1);
   }

   // the property definitions of an object or of a context (the project)
   std::vector<IfcSchema::IfcRelDefinesByProperties> defined_by(IfcSchema::IfcObjectDefinition object)
   {
      if (auto o = object.as<IfcSchema::IfcObject>())
         return o.IsDefinedBy();
      if (auto c = object.as<IfcSchema::IfcContext>())
         return c.IsDefinedBy();
      return {};
   }

   // The world transformation of a placement (lengths in meters), or nullopt if it can't be computed here (e.g. a linear
   // placement without a Cartesian position)
   std::optional<Eigen::Matrix4d> axis2placement_matrix(const express::base& placement, double meters)
   {
      Eigen::Matrix4d m = Eigen::Matrix4d::Identity();
      if (auto a3 = placement.as<IfcSchema::IfcAxis2Placement3D>())
      {
         auto point = a3.Location().as<IfcSchema::IfcCartesianPoint>();
         if (!point)
            return std::nullopt;
         auto c = point.Coordinates();
         Eigen::Vector3d origin(c.size() > 0 ? c[0] : 0, c.size() > 1 ? c[1] : 0, c.size() > 2 ? c[2] : 0);
         Eigen::Vector3d z(0, 0, 1), x(1, 0, 0);
         if (auto axis = a3.Axis())
         {
            auto r = axis.DirectionRatios();
            z = Eigen::Vector3d(r[0], r[1], r.size() > 2 ? r[2] : 0).normalized();
         }
         if (auto ref = a3.RefDirection())
         {
            auto r = ref.DirectionRatios();
            x = Eigen::Vector3d(r[0], r[1], r.size() > 2 ? r[2] : 0);
         }
         else if (IsEqual(fabs(z.z()), 1.0))
         {
            x = Eigen::Vector3d(1, 0, 0);
         }
         else
         {
            x = Eigen::Vector3d(0, 0, 1).cross(z);
         }
         x = (x - x.dot(z) * z).normalized();
         Eigen::Vector3d y = z.cross(x);
         m.block<3, 1>(0, 0) = x;
         m.block<3, 1>(0, 1) = y;
         m.block<3, 1>(0, 2) = z;
         m.block<3, 1>(0, 3) = origin * meters;
         return m;
      }
      if (auto a2 = placement.as<IfcSchema::IfcAxis2Placement2D>())
      {
         auto point = a2.Location().as<IfcSchema::IfcCartesianPoint>();
         if (!point)
            return std::nullopt;
         auto c = point.Coordinates();
         Eigen::Vector3d x(1, 0, 0);
         if (auto ref = a2.RefDirection())
         {
            auto r = ref.DirectionRatios();
            x = Eigen::Vector3d(r[0], r[1], 0).normalized();
         }
         m.block<3, 1>(0, 0) = x;
         m.block<3, 1>(0, 1) = Eigen::Vector3d(-x.y(), x.x(), 0);
         m.block<3, 1>(0, 3) = Eigen::Vector3d(c[0], c.size() > 1 ? c[1] : 0, 0) * meters;
         return m;
      }
      return std::nullopt;
   }

   std::optional<Eigen::Matrix4d> placement_matrix(IfcSchema::IfcObjectPlacement placement, double meters, int depth = 0)
   {
      if (!placement)
         return Eigen::Matrix4d::Identity();
      if (100 < depth)
         return std::nullopt; // a cycle

      std::optional<Eigen::Matrix4d> relative;
      if (auto local = placement.as<IfcSchema::IfcLocalPlacement>())
         relative = axis2placement_matrix(local.RelativePlacement(), meters);
      else if (auto linear = placement.as<IfcSchema::IfcLinearPlacement>(); linear && linear.CartesianPosition())
         relative = axis2placement_matrix(linear.CartesianPosition(), meters);
      if (!relative)
         return std::nullopt;

      auto parent = placement_matrix(placement.PlacementRelTo(), meters, depth + 1);
      if (!parent)
         return std::nullopt;
      return (*parent) * (*relative);
   }

   // The relationships the merge rebuilds for the elements it adds: the relating and related attributes
   struct RelationshipKind
   {
      const char* type;
      const char* relating;
      const char* related;
   };
   const std::array<RelationshipKind, 7> relationship_kinds{ {
      {"IfcRelAggregates", "RelatingObject", "RelatedObjects"},
      {"IfcRelNests", "RelatingObject", "RelatedObjects"},
      {"IfcRelContainedInSpatialStructure", "RelatingStructure", "RelatedElements"},
      {"IfcRelDefinesByType", "RelatingType", "RelatedObjects"},
      {"IfcRelDefinesByProperties", "RelatingPropertyDefinition", "RelatedObjects"},
      {"IfcRelAssociatesMaterial", "RelatingMaterial", "RelatedObjects"},
      {"IfcRelAssociatesClassification", "RelatingClassification", "RelatedObjects"},
   } };

   class Merger
   {
   public:
      Merger(ifcopenshell::file& exported, ifcopenshell::file& source, const IfcExportIdentities& identities, const CRoundTripRegistry& registry, const CIfcMappingTable& table, RoundTripMergeResult& result) :
         m_Exp(exported), m_Src(source), m_Identities(identities), m_Registry(registry), m_Table(table), m_Result(result)
      {
         m_ExpUnits.Init(m_Exp);
         m_SrcUnits.Init(m_Src);
         m_ExpMeters = m_ExpUnits.GetConversionFactor({}, IfcSchema::IfcUnitEnum::IfcUnit_LENGTHUNIT);
         m_SrcMeters = m_SrcUnits.GetConversionFactor({}, IfcSchema::IfcUnitEnum::IfcUnit_LENGTHUNIT);
         m_SourceMaxId = m_Src.get_max_id();
         for (const auto& context : m_Src.instances_by_type<IfcSchema::IfcRepresentationContext>())
            m_SourceContexts.insert(context.id());
         for (const auto& classification : m_Src.instances_by_type<IfcSchema::IfcClassification>())
            m_SourceClassifications.insert(classification.id());
      }

      void Run();

   private:
      ifcopenshell::file& m_Exp;
      ifcopenshell::file& m_Src;
      const IfcExportIdentities& m_Identities;
      const CRoundTripRegistry& m_Registry;
      const CIfcMappingTable& m_Table;
      RoundTripMergeResult& m_Result;
      CIfcImportUnits m_ExpUnits, m_SrcUnits;
      double m_ExpMeters = 1.0, m_SrcMeters = 1.0;
      unsigned int m_SourceMaxId = 0;
      std::set<unsigned int> m_SourceContexts, m_SourceClassifications;

      struct Match
      {
         std::string key;
         IfcSchema::IfcObjectDefinition exported, source;
      };
      std::vector<Match> m_Matches;
      std::map<unsigned int, IfcSchema::IfcObjectDefinition> m_SourceOf; // exported element id -> source element

      // an element of PGSuper's model to add, and the source element its placement is relative to (may be empty).
      // When into is set, only the element's parts are added, to that source element (e.g. strands to an empty bundle)
      struct Addition
      {
         IfcSchema::IfcProduct exported;
         IfcSchema::IfcProduct anchor;
         std::string description;
         IfcSchema::IfcObjectDefinition into;
      };
      std::vector<Addition> m_Additions;

      // the values given to shared type and material properties, so a different value makes a copy (see Assign)
      std::map<std::pair<unsigned int, std::string>, std::string> m_SharedValues;
      std::map<std::pair<unsigned int, std::string>, IfcSchema::IfcMaterialDefinition> m_MaterialCopies; // the copies of a source material, by the values given to them

      void Log(const std::string& line) { m_Result.log.push_back(line); }
      static std::string Describe(const express::base& element);

      IfcSchema::IfcObjectDefinition FindSource(const std::string& global_id);
      void MatchElements();
      void PlanChildren(const Match& girder);
      ElementKind ChildRole(const IfcSchema::IfcObject& child, const CIfcTargetReader& reader) const;

      void AddElements();
      void FixCopiedValues(const std::set<unsigned int>& copied);
      void Redirect(const express::base& from, const express::base& to);
      void RemoveIfUnreferenced(const express::base& entity, int depth = 0);
      bool IsNew(const express::base& entity) const { return m_SourceMaxId < entity.id(); }
      void MergeCopiedContextsAndClassifications();

      void UpdateElement(const Match& match);
      void UpdateClassifications(const Match& match);
      const PropertyDeclaration* FindDeclaration(const std::string& key, PropertyOwner attach, const std::string& pset, const std::string& property) const;

      // Updates the source's property set (a property set, quantity set, or material properties) from PGSuper's. Returns the
      // number of properties changed and added
      std::pair<IndexType, IndexType> UpdateProperties(const std::string& key, PropertyOwner attach, const express::base& exported_pset, express::entity source_pset, const std::string& where, bool bShared);
      express::base ConvertedValue(const express::base& value, IfcSchema::IfcUnit from_unit, IfcSchema::IfcUnit to_unit);
      express::base CopyProperty(const express::base& property);
      express::base CopyQuantity(const express::base& quantity);

      IfcSchema::IfcPropertySetDefinition FindOccurrencePset(IfcSchema::IfcObjectDefinition object, const std::string& name, IfcSchema::IfcRelDefinesByProperties* pRel = nullptr);
      IfcSchema::IfcPropertySetDefinition OwnOccurrencePset(IfcSchema::IfcObjectDefinition object, IfcSchema::IfcPropertySetDefinition pset, IfcSchema::IfcRelDefinesByProperties rel);
      IfcSchema::IfcPropertySetDefinition NewOccurrencePset(IfcSchema::IfcObjectDefinition object, const express::base& exported_pset);
   };

   std::string Merger::Describe(const express::base& element)
   {
      std::ostringstream os;
      os << element.declaration().name();
      if (auto root = element.as<IfcSchema::IfcRoot>())
         os << " \"" << root.Name().value_or("unnamed") << "\" (" << root.GlobalId() << ")";
      else
         os << " #" << element.id();
      return os.str();
   }

   IfcSchema::IfcObjectDefinition Merger::FindSource(const std::string& global_id)
   {
      try
      {
         return m_Src.instance_by_guid(global_id).as<IfcSchema::IfcObjectDefinition>();
      }
      catch (...)
      {
         return {};
      }
   }

   void Merger::MatchElements()
   {
      for (const auto& [key, exported_ids] : m_Identities)
      {
         std::vector<IfcSchema::IfcObjectDefinition> exported;
         for (auto id : exported_ids)
         {
            if (auto e = m_Exp.instance_by_id(id).as<IfcSchema::IfcObjectDefinition>())
               exported.push_back(e);
         }
         std::vector<IfcSchema::IfcObjectDefinition> source;
         for (const auto& global_id : m_Registry.Find(key))
         {
            if (auto s = FindSource(global_id))
               source.push_back(s);
         }

         auto prefix = key.substr(0, key.find(':'));
         if (prefix == "bearings")
         {
            if (source.empty())
            {
               for (const auto& e : exported)
                  m_Additions.push_back({ e.as<IfcSchema::IfcProduct>(), {}, "bearing" });
            }
            else if (source.size() == exported.size())
            {
               // pair each of PGSuper's bearings with the nearest of the source's (by placement origin), else in order
               std::vector<bool> used(source.size(), false);
               for (const auto& e : exported)
               {
                  auto me = placement_matrix(e.as<IfcSchema::IfcProduct>().ObjectPlacement(), m_ExpMeters);
                  size_t best = source.size();
                  double best_distance = DBL_MAX;
                  for (size_t i = 0; i < source.size(); i++)
                  {
                     if (used[i])
                        continue;
                     auto ms = placement_matrix(source[i].as<IfcSchema::IfcProduct>().ObjectPlacement(), m_SrcMeters);
                     double d = (me && ms) ? ((*me).block<3, 1>(0, 3) - (*ms).block<3, 1>(0, 3)).norm() : (double)i;
                     if (d < best_distance)
                     {
                        best_distance = d;
                        best = i;
                     }
                  }
                  used[best] = true;
                  m_Matches.push_back({ key, e, source[best] });
               }
            }
            else
            {
               std::ostringstream os;
               os << "The model has " << source.size() << " bearings at " << key << " and PGSuper " << exported.size() << ": the model's bearings aren't changed.";
               Log(os.str());
            }
            continue;
         }

         if (prefix == "diaphragms")
         {
            if (source.empty())
            {
               for (const auto& e : exported)
                  m_Additions.push_back({ e.as<IfcSchema::IfcProduct>(), {}, "pier diaphragm" });
            }
            else if (!exported.empty())
            {
               // the model's diaphragms at the pier get PGSuper's diaphragm values; PGSuper's aren't added
               for (const auto& s : source)
                  m_Matches.push_back({ key, exported.front(), s });
            }
            continue;
         }

         if (exported.empty())
            continue;

         if (source.empty())
         {
            if (prefix == "girder" || prefix == "deck" || prefix == "barrier" || prefix == "pier")
               Log(Describe(exported.front()) + " of PGSuper's model isn't in the IFC model, so it isn't added (the round trip doesn't add girders, decks, barriers, or piers yet).");
            continue;
         }

         m_Matches.push_back({ key, exported.front(), source.front() });
      }

      for (const auto& match : m_Matches)
      {
         m_SourceOf[match.exported.id()] = match.source;
         if (auto r1 = match.exported.as<IfcSchema::IfcRoot>(); r1)
            m_Result.global_ids[r1.GlobalId()] = match.source.as<IfcSchema::IfcRoot>().GlobalId();
      }

      // PlanChildren can add matches, so the girder is a copy
      for (size_t i = 0, n = m_Matches.size(); i < n; i++)
      {
         if (m_Matches[i].key.rfind("girder:", 0) == 0)
         {
            Match girder = m_Matches[i];
            PlanChildren(girder);
         }
      }
   }

   // The element role of a child of a girder, as an assembly of strands or bars: by its parts, else by the mapping table
   // (e.g. an empty assembly). ElementKind::Girder if it's neither
   ElementKind Merger::ChildRole(const IfcSchema::IfcObject& child, const CIfcTargetReader& reader) const
   {
      bool bTendons = false, bBars = false;
      for (const auto& part_rel : child.IsDecomposedBy())
      {
         for (const auto& part : part_rel.RelatedObjects())
         {
            bTendons |= (bool)part.as<IfcSchema::IfcTendon>();
            bBars |= (bool)part.as<IfcSchema::IfcReinforcingBar>();
         }
      }
      if (bTendons != bBars)
         return bTendons ? ElementKind::TendonBundle : ElementKind::ReinforcementCage;
      if (bTendons)
         return ElementKind::Girder; // both: not one of PGSuper's

      if (reader.Matches(ElementKind::ReinforcementCage, child))
         return ElementKind::ReinforcementCage;
      if (reader.Matches(ElementKind::TendonBundle, child))
         return ElementKind::TendonBundle;
      return ElementKind::Girder;
   }

   // PGSuper's tendon bundles and reinforcement cages of a girder are matched with the source girder's of the same role,
   // which are then updated like other matched elements (property sets, classification). PGSuper's strands and bars are
   // added if the source's girder has none: into the matched assembly if it's empty, else with their own assembly
   void Merger::PlanChildren(const Match& girder)
   {
      auto segment = id_of_key(girder.key);
      bool bSourceTendons = !m_Registry.Find("tendons:" + segment).empty();
      bool bSourceBars = !m_Registry.Find("bars:" + segment).empty();
      std::string name = girder.source.as<IfcSchema::IfcRoot>().Name().value_or("unnamed girder");

      CIfcTargetReader exported_reader(m_Table, m_ExpUnits), source_reader(m_Table, m_SrcUnits);
      auto has_parts = [](const IfcSchema::IfcObjectDefinition& object)
         {
            for (const auto& rel : object.IsDecomposedBy())
               if (!rel.RelatedObjects().empty())
                  return true;
            for (const auto& rel : object.IsNestedBy())
               if (!rel.RelatedObjects().empty())
                  return true;
            return false;
         };

      std::set<unsigned int> used; // the source girder's children already matched
      for (const auto& rel : girder.exported.IsDecomposedBy())
      {
         for (const auto& child : rel.RelatedObjects())
         {
            auto product = child.as<IfcSchema::IfcProduct>();
            if (!product)
               continue;

            auto role = ChildRole(product, exported_reader);
            if (role == ElementKind::Girder)
               continue;
            bool bTendons = (role == ElementKind::TendonBundle);
            bool bSourceParts = (bTendons ? bSourceTendons : bSourceBars);
            std::string key = (bTendons ? "tendons:" : "bars:") + segment;

            // the source girder's assembly of the role: an empty one, or the one with the model's strands or bars
            IfcSchema::IfcObjectDefinition source_child;
            for (const auto& source_rel : girder.source.IsDecomposedBy())
            {
               for (const auto& candidate : source_rel.RelatedObjects())
               {
                  auto object = candidate.as<IfcSchema::IfcObject>();
                  if (source_child || !object || used.count(object.id()) || ChildRole(object, source_reader) != role)
                     continue;
                  if (!has_parts(object) || bSourceParts)
                     source_child = object;
               }
            }

            if (source_child)
            {
               used.insert(source_child.id());
               m_Matches.push_back({ key, product, source_child });
               m_SourceOf[product.id()] = source_child;
               m_Result.global_ids[product.GlobalId()] = source_child.as<IfcSchema::IfcRoot>().GlobalId();
            }

            if (bSourceParts || !has_parts(product))
               continue;

            // into the source's assembly only if it's empty (one with other parts isn't used)
            m_Additions.push_back({ product, girder.source.as<IfcSchema::IfcProduct>(), std::string(bTendons ? "strands" : "bars") + " of " + name, source_child });
         }
      }

      if (bSourceTendons)
         Log("Girder \"" + name + "\" keeps the model's strands; PGSuper's aren't added.");
      if (bSourceBars)
         Log("Girder \"" + name + "\" keeps the model's reinforcing bars; PGSuper's aren't added.");
      m_Result.bModelTendons |= bSourceTendons;
      m_Result.bModelBars |= bSourceBars;
   }

   void Merger::Redirect(const express::base& from, const express::base& to)
   {
      for (auto& instance : m_Src.instances_by_reference(from.id()))
      {
         auto nAttributes = instance.declaration().as_entity() ? instance.declaration().as_entity()->attribute_count() : 0;
         for (size_t i = 0; i < nAttributes; i++)
         {
            auto attribute = instance.get_attribute_value(i);
            switch (attribute.type())
            {
            case ifcopenshell::Argument_ENTITY_INSTANCE:
               if ((express::base)attribute == from)
                  instance.set_attribute_value(i, to);
               break;
            case ifcopenshell::Argument_AGGREGATE_OF_ENTITY_INSTANCE:
            {
               std::vector<express::base> list = attribute;
               if (std::find(list.begin(), list.end(), from) != list.end())
               {
                  std::replace(list.begin(), list.end(), from, to);
                  instance.set_attribute_value(i, list);
               }
               break;
            }
            default:
               break;
            }
         }
      }
   }

   // Removes an entity the merge added if nothing refers to it any more, and then what it refers to (never the source's own)
   void Merger::RemoveIfUnreferenced(const express::base& entity, int depth)
   {
      if (!entity || !IsNew(entity) || 0 < m_Src.get_total_inverses(entity.id()) || 200 < depth)
         return;

      std::vector<express::base> references;
      auto nAttributes = entity.declaration().as_entity() ? entity.declaration().as_entity()->attribute_count() : 0;
      for (size_t i = 0; i < nAttributes; i++)
      {
         auto attribute = entity.get_attribute_value(i);
         if (attribute.type() == ifcopenshell::Argument_ENTITY_INSTANCE)
            references.push_back(attribute);
         else if (attribute.type() == ifcopenshell::Argument_AGGREGATE_OF_ENTITY_INSTANCE)
         {
            std::vector<express::base> list = attribute;
            references.insert(references.end(), list.begin(), list.end());
         }
      }
      m_Src.remove_entity(entity);
      for (const auto& reference : references)
      {
         if (reference.declaration().as_entity())
            RemoveIfUnreferenced(reference, depth + 1);
      }
   }

   void Merger::AddElements()
   {
      if (m_Additions.empty())
         return;

      // the elements to add, with their parts, and the source element each one's placement is relative to
      std::set<unsigned int> elements;
      std::map<unsigned int, IfcSchema::IfcProduct> anchor_of;
      std::vector<std::pair<IfcSchema::IfcObjectDefinition, IfcSchema::IfcProduct>> stack;
      for (const auto& addition : m_Additions)
      {
         if (!addition.into)
         {
            stack.emplace_back(addition.exported, addition.anchor);
            continue;
         }
         // only the parts: the source element takes the place of PGSuper's (see m_SourceOf)
         for (const auto& rel : addition.exported.IsDecomposedBy())
            for (const auto& part : rel.RelatedObjects())
               stack.emplace_back(part, addition.anchor);
         for (const auto& rel : addition.exported.IsNestedBy())
            for (const auto& part : rel.RelatedObjects())
               stack.emplace_back(part, addition.anchor);
      }
      while (!stack.empty())
      {
         auto [element, anchor] = stack.back();
         stack.pop_back();
         if (!elements.insert(element.id()).second)
            continue;
         anchor_of[element.id()] = anchor;
         for (const auto& rel : element.IsDecomposedBy())
            for (const auto& part : rel.RelatedObjects())
               stack.emplace_back(part, anchor);
         for (const auto& rel : element.IsNestedBy())
            for (const auto& part : rel.RelatedObjects())
               stack.emplace_back(part, anchor);
      }

      // copy them (IfcOpenShell converts lengths to the source's unit), and note everything copied with them
      std::set<unsigned int> copied;
      std::map<unsigned int, express::base> copy_of;
      auto copy = [&](const express::base& e)
         {
            for (const auto& reached : ifcopenshell::file::traverse(e))
            {
               if (reached.id() != 0) // values (labels, measures) are inline: id 0
                  copied.insert(reached.id());
            }
            auto c = m_Src.add_entity(e);
            copy_of[e.id()] = c;
            return c;
         };
      for (auto id : elements)
         copy(m_Exp.instance_by_id(id));

      // Place the added elements in the source. A placement relative to another added element's stays as it is; any other
      // (relative to an element of PGSuper's model the source doesn't share, e.g. PGSuper's girder, or a linear placement
      // along PGSuper's alignment) is made relative to the source element the added element belongs to, when its placement
      // can be computed, else to the world
      std::set<unsigned int> own_placements;
      for (auto id : elements)
      {
         if (auto product = m_Exp.instance_by_id(id).as<IfcSchema::IfcProduct>(); product && product.ObjectPlacement())
            own_placements.insert(product.ObjectPlacement().id());
      }
      std::set<unsigned int> placed;
      std::vector<express::base> replaced; // what the new placements replaced, removed when nothing refers to it
      std::set<unsigned int> unplaced;
      for (auto id : elements)
      {
         auto product = m_Exp.instance_by_id(id).as<IfcSchema::IfcProduct>();
         if (!product || !product.ObjectPlacement())
            continue;
         auto placement = product.ObjectPlacement();
         auto parent = placement.PlacementRelTo();
         bool bLocal = (bool)placement.as<IfcSchema::IfcLocalPlacement>();
         if ((bLocal && parent && own_placements.count(parent.id())) || !placed.insert(placement.id()).second)
            continue;

         auto world = placement_matrix(placement, m_ExpMeters);
         if (!world)
         {
            unplaced.insert(id);
            continue;
         }

         IfcSchema::IfcObjectPlacement relative_to;
         Eigen::Matrix4d m = *world;
         auto anchor = anchor_of[id];
         if (anchor && anchor.ObjectPlacement())
         {
            if (auto anchor_matrix = placement_matrix(anchor.ObjectPlacement(), m_SrcMeters))
            {
               relative_to = anchor.ObjectPlacement();
               m = anchor_matrix->inverse() * (*world);
            }
         }

         auto location = m_Src.create<IfcSchema::IfcCartesianPoint>().initialize(std::vector<double>{ m(0, 3) / m_SrcMeters, m(1, 3) / m_SrcMeters, m(2, 3) / m_SrcMeters });
         auto axis = m_Src.create<IfcSchema::IfcDirection>().initialize(std::vector<double>{ m(0, 2), m(1, 2), m(2, 2) });
         auto ref = m_Src.create<IfcSchema::IfcDirection>().initialize(std::vector<double>{ m(0, 0), m(1, 0), m(2, 0) });
         auto a2p = m_Src.create<IfcSchema::IfcAxis2Placement3D>().initialize(location, axis, ref);

         auto copied_placement = copy_of[id].as<IfcSchema::IfcProduct>().ObjectPlacement();
         if (auto local = copied_placement.as<IfcSchema::IfcLocalPlacement>())
         {
            // changed in place, so the placements relative to it follow
            replaced.push_back(local.RelativePlacement());
            replaced.push_back(local.PlacementRelTo());
            local.setPlacementRelTo(relative_to);
            local.setRelativePlacement(a2p);
         }
         else
         {
            auto new_placement = m_Src.create<IfcSchema::IfcLocalPlacement>().initialize(relative_to, a2p);
            Redirect(copied_placement, new_placement);
            replaced.push_back(copied_placement);
         }
      }
      for (const auto& entity : replaced)
         RemoveIfUnreferenced(entity);
      if (!unplaced.empty())
         Log(std::to_string(unplaced.size()) + " added elements have placements that can't be computed here; they keep PGSuper's placements.");

      for (const auto& addition : m_Additions)
      {
         m_Result.nAdded++;
         if (addition.into)
            Log("Added " + addition.description + " to the model's empty " + Describe(addition.into) + ".");
         else
            Log("Added " + addition.description + ": " + Describe(addition.exported) + ".");
      }

      // the relationships of the added elements, with only the added elements in them. Types, materials, property sets, and
      // classifications are copied; an element of PGSuper's model is replaced by its source element
      std::set<unsigned int> related_set = elements; // the added elements, then also their types
      auto map_related = [&](const express::base& e) -> express::base
         {
            auto found = copy_of.find(e.id());
            return found == copy_of.end() ? express::base() : found->second;
         };
      auto rebuild = [&](const RelationshipKind& kind, bool bTypes)
         {
            for (const auto& rel : m_Exp.instances_by_type(kind.type))
            {
               auto entity = rel.as<express::entity>();
               std::vector<express::base> related = entity.get(kind.related);
               std::vector<express::base> new_related;
               for (const auto& r : related)
               {
                  if (related_set.count(r.id()))
                  {
                     auto c = map_related(r);
                     if (!c)
                        c = copy(r);
                     new_related.push_back(c);
                  }
               }
               if (new_related.empty())
                  continue;

               express::base relating = entity.get(kind.relating);
               express::base new_relating;
               if (auto found = copy_of.find(relating.id()); found != copy_of.end())
                  new_relating = found->second;
               else if (auto source = m_SourceOf.find(relating.id()); source != m_SourceOf.end())
                  new_relating = source->second;
               else if (relating.as<IfcSchema::IfcProduct>() || relating.as<IfcSchema::IfcSpatialElement>())
               {
                  // a spatial container the source doesn't have: the source's bridge
                  auto bridge = m_Registry.Find(roundtrip_key::Bridge());
                  if (!bridge.empty() && std::string(kind.type) == "IfcRelContainedInSpatialStructure")
                     new_relating = FindSource(bridge.front());
                  if (!new_relating)
                  {
                     Log("A relationship of the added elements to " + Describe(relating) + " isn't in the IFC model, so it's left out.");
                     continue;
                  }
               }
               else
               {
                  new_relating = copy(relating);
               }

               if (bTypes)
                  related_set.insert(relating.id()); // the type's materials and classifications are related to it

               auto new_rel = m_Src.create(&rel.declaration());
               new_rel.set_attribute_value("GlobalId", (std::string)ifcopenshell::global_id());
               if (auto name = entity.get("Name"); !name.isNull())
                  new_rel.set_attribute_value("Name", (std::string)name);
               new_rel.set_attribute_value(kind.relating, new_relating);
               new_rel.set_attribute_value(kind.related, new_related);
            }
         };
      for (const auto& kind : relationship_kinds)
      {
         std::string type(kind.type);
         if (type == "IfcRelAggregates" || type == "IfcRelNests" || type == "IfcRelContainedInSpatialStructure")
            rebuild(kind, false);
      }
      for (const auto& kind : relationship_kinds)
      {
         if (std::string(kind.type) == "IfcRelDefinesByType")
            rebuild(kind, true);
      }
      for (const auto& kind : relationship_kinds)
      {
         std::string type(kind.type);
         if (type == "IfcRelDefinesByProperties" || type == "IfcRelAssociatesMaterial" || type == "IfcRelAssociatesClassification")
            rebuild(kind, false);
      }

      // what refers to the copies: the styles of their geometry and the properties of their materials (the export's materials
      // have no representations)
      for (const auto& styled_item : m_Exp.instances_by_type<IfcSchema::IfcStyledItem>())
      {
         if (styled_item.Item() && copied.count(styled_item.Item().id()))
            copy(styled_item);
      }
      for (const auto& properties : m_Exp.instances_by_type<IfcSchema::IfcMaterialProperties>())
      {
         if (properties.Material() && copied.count(properties.Material().id()))
            copy(properties);
      }

      FixCopiedValues(copied);
      MergeCopiedContextsAndClassifications();
   }

   // IfcOpenShell converts the lengths of copied attributes typed IfcLengthMeasure only. Property and quantity values are
   // converted here, from PGSuper's units to the source's (a value with its own unit keeps it)
   void Merger::FixCopiedValues(const std::set<unsigned int>& copied)
   {
      for (auto id : copied)
      {
         auto exported = m_Exp.instance_by_id(id);
         if (!exported)
            continue;
         if (auto property = exported.as<IfcSchema::IfcPropertySingleValue>())
         {
            if (!property.NominalValue() || property.Unit())
               continue;
            auto c = m_Src.add_entity(exported).as<IfcSchema::IfcPropertySingleValue>(); // the copy made before
            if (auto value = ConvertedValue(property.NominalValue(), {}, {}))
               c.setNominalValue(IfcSchema::IfcValue(value));
         }
         else if (exported.as<IfcSchema::IfcPhysicalSimpleQuantity>())
         {
            CopyQuantity(exported); // the copy made before, its value converted
         }
      }
   }

   void Merger::MergeCopiedContextsAndClassifications()
   {
      // representation contexts: the source's, by identifier and type
      std::vector<IfcSchema::IfcRepresentationContext> source_contexts;
      for (auto id : m_SourceContexts)
         source_contexts.push_back(m_Src.instance_by_id(id).as<IfcSchema::IfcRepresentationContext>());
      auto main_context = [&]() -> IfcSchema::IfcRepresentationContext
         {
            for (const auto& c : source_contexts)
            {
               if (c.as<IfcSchema::IfcGeometricRepresentationContext>() && !c.as<IfcSchema::IfcGeometricRepresentationSubContext>() && c.ContextType().value_or("") == "Model")
                  return c;
            }
            for (const auto& c : source_contexts)
            {
               if (c.as<IfcSchema::IfcGeometricRepresentationContext>() && !c.as<IfcSchema::IfcGeometricRepresentationSubContext>())
                  return c;
            }
            return {};
         };

      for (const auto& context : m_Src.instances_by_type<IfcSchema::IfcRepresentationContext>())
      {
         if (m_SourceContexts.count(context.id()))
            continue;

         IfcSchema::IfcRepresentationContext replacement;
         for (const auto& c : source_contexts)
         {
            if (c.ContextIdentifier() == context.ContextIdentifier() && c.ContextType() == context.ContextType())
            {
               replacement = c;
               break;
            }
         }
         if (!replacement && context.as<IfcSchema::IfcGeometricRepresentationSubContext>())
         {
            // the source has no such sub-context: make one in its main context
            if (auto parent = main_context())
            {
               auto sub = context.as<IfcSchema::IfcGeometricRepresentationSubContext>();
               auto new_sub = m_Src.create<IfcSchema::IfcGeometricRepresentationSubContext>();
               new_sub.setContextIdentifier(sub.ContextIdentifier());
               new_sub.setContextType(sub.ContextType());
               new_sub.setParentContext(parent.as<IfcSchema::IfcGeometricRepresentationContext>());
               new_sub.setTargetView(sub.TargetView());
               replacement = new_sub;
               source_contexts.push_back(new_sub);
               m_SourceContexts.insert(new_sub.id());
            }
         }
         if (!replacement)
            replacement = main_context();
         if (replacement)
            Redirect(context, replacement);
      }
      for (const auto& context : m_Src.instances_by_type<IfcSchema::IfcRepresentationContext>())
      {
         if (!m_SourceContexts.count(context.id()))
            RemoveIfUnreferenced(context);
      }

      // classification systems: the source's of the same name, else the copy, associated with the source's project
      for (const auto& classification : m_Src.instances_by_type<IfcSchema::IfcClassification>())
      {
         if (m_SourceClassifications.count(classification.id()))
            continue;

         IfcSchema::IfcClassification existing;
         for (auto id : m_SourceClassifications)
         {
            auto c = m_Src.instance_by_id(id).as<IfcSchema::IfcClassification>();
            if (c && c.Name() == classification.Name())
            {
               existing = c;
               break;
            }
         }
         if (existing)
         {
            Redirect(classification, existing);
            RemoveIfUnreferenced(classification);
         }
         else
         {
            m_SourceClassifications.insert(classification.id());
            auto projects = m_Src.instances_by_type<IfcSchema::IfcProject>();
            if (!projects.empty())
            {
               auto rel = m_Src.create<IfcSchema::IfcRelAssociatesClassification>();
               rel.setGlobalId(ifcopenshell::global_id());
               rel.setRelatedObjects(std::vector<IfcSchema::IfcDefinitionSelect>{ projects.front() });
               rel.setRelatingClassification(classification);
            }
            Log("Added the classification system \"" + classification.Name() + "\".");
         }
      }
   }

   const PropertyDeclaration* Merger::FindDeclaration(const std::string& key, PropertyOwner attach, const std::string& pset, const std::string& property) const
   {
      for (auto role : roles_of_key(key))
      {
         for (const auto* declaration : m_Table.GetPropertySets(role, attach))
         {
            if (declaration->name != pset)
               continue;
            for (const auto& p : declaration->properties)
            {
               if (p.name == property)
                  return &p;
            }
         }
      }
      return nullptr;
   }

   // A value in the source's units: from from_unit (else PGSuper's project unit) to to_unit (else the source's project unit).
   // Values without a unit are copied
   express::base Merger::ConvertedValue(const express::base& value, IfcSchema::IfcUnit from_unit, IfcSchema::IfcUnit to_unit)
   {
      if (!value)
         return {};
      auto unit_type = unit_type_of_measure(value.declaration().name());
      auto number = number_of(value);
      if (!unit_type || !number)
         return m_Src.add_entity(value);

      double from = m_ExpUnits.GetConversionFactor(from_unit, *unit_type);
      double to = m_SrcUnits.GetConversionFactor(to_unit, *unit_type);
      auto converted = m_Src.create(&value.declaration());
      converted.set_attribute_value(0, *number * from / to);
      return converted;
   }

   // A copy of PGSuper's property for the source, its value in the source's project unit (a value with its own unit keeps it)
   express::base Merger::CopyProperty(const express::base& property)
   {
      auto single = property.as<IfcSchema::IfcPropertySingleValue>();
      if (!single || single.Unit() || !single.NominalValue())
         return m_Src.add_entity(property);

      auto copy = m_Src.create<IfcSchema::IfcPropertySingleValue>();
      copy.setName(single.Name());
      copy.setSpecification(single.Specification());
      copy.setNominalValue(IfcSchema::IfcValue(ConvertedValue(single.NominalValue(), {}, {})));
      return copy;
   }

   // A copy of PGSuper's quantity for the source, its value in the source's project unit (a value with its own unit keeps it)
   express::base Merger::CopyQuantity(const express::base& quantity)
   {
      auto copy = m_Src.add_entity(quantity); // converts a length; the other quantities are converted here
      auto value_attribute = quantity_value_of(quantity.declaration().name());
      if (!value_attribute)
         return copy;
      auto unit = IfcSchema::IfcUnit(quantity.as<IfcSchema::IfcPhysicalSimpleQuantity>().Unit());
      double v = quantity.as<express::entity>().get(value_attribute->first);
      double from = m_ExpUnits.GetConversionFactor(unit, value_attribute->second);
      double to = unit ? from : m_SrcUnits.GetConversionFactor({}, value_attribute->second);
      copy.set_attribute_value(value_attribute->first, v * from / to);
      return copy;
   }

   IfcSchema::IfcPropertySetDefinition Merger::FindOccurrencePset(IfcSchema::IfcObjectDefinition object, const std::string& name, IfcSchema::IfcRelDefinesByProperties* pRel)
   {
      for (const auto& rel : defined_by(object))
      {
         auto definition = rel.RelatingPropertyDefinition().as<IfcSchema::IfcPropertySetDefinition>();
         if (definition && definition.Name().value_or("") == name)
         {
            if (pRel)
               *pRel = rel;
            return definition;
         }
      }
      return {};
   }

   // The property set as the element's own: a property set related to other elements too is copied for the element
   IfcSchema::IfcPropertySetDefinition Merger::OwnOccurrencePset(IfcSchema::IfcObjectDefinition object, IfcSchema::IfcPropertySetDefinition pset, IfcSchema::IfcRelDefinesByProperties rel)
   {
      size_t nObjects = 0;
      for (const auto& r : pset.DefinesOccurrence())
         nObjects += r.RelatedObjects().size();
      if (nObjects <= 1)
         return pset;

      // a copy of the property set with the same properties. A property that changes is replaced in the copy only, since
      // it's shared (see UpdateProperties)
      auto copy = m_Src.create(&pset.declaration());
      auto entity = pset.as<express::entity>();
      copy.set_attribute_value("GlobalId", (std::string)ifcopenshell::global_id());
      if (auto name = entity.get("Name"); !name.isNull())
         copy.set_attribute_value("Name", (std::string)name);
      if (auto description = entity.get("Description"); !description.isNull())
         copy.set_attribute_value("Description", (std::string)description);
      if (pset.as<IfcSchema::IfcElementQuantity>())
      {
         if (auto method = entity.get("MethodOfMeasurement"); !method.isNull())
            copy.set_attribute_value("MethodOfMeasurement", (std::string)method);
         copy.set_attribute_value("Quantities", (std::vector<express::base>)entity.get("Quantities"));
      }
      else
      {
         copy.set_attribute_value("HasProperties", (std::vector<express::base>)entity.get("HasProperties"));
      }

      // the element leaves the shared relationship and gets its own
      auto related = rel.RelatedObjects();
      related.erase(std::remove(related.begin(), related.end(), object), related.end());
      rel.setRelatedObjects(related);
      auto own = m_Src.create<IfcSchema::IfcRelDefinesByProperties>();
      own.setGlobalId(ifcopenshell::global_id());
      own.setRelatedObjects(std::vector<IfcSchema::IfcObjectDefinition>{ object });
      own.setRelatingPropertyDefinition(IfcSchema::IfcPropertySetDefinitionSelect(copy));
      return copy.as<IfcSchema::IfcPropertySetDefinition>();
   }

   // A new property set (or quantity set) on the element, empty, with the name and kind of PGSuper's
   IfcSchema::IfcPropertySetDefinition Merger::NewOccurrencePset(IfcSchema::IfcObjectDefinition object, const express::base& exported_pset)
   {
      auto pset = m_Src.create(&exported_pset.declaration());
      auto entity = exported_pset.as<express::entity>();
      pset.set_attribute_value("GlobalId", (std::string)ifcopenshell::global_id());
      if (auto name = entity.get("Name"); !name.isNull())
         pset.set_attribute_value("Name", (std::string)name);
      if (auto description = entity.get("Description"); !description.isNull())
         pset.set_attribute_value("Description", (std::string)description);
      if (exported_pset.as<IfcSchema::IfcElementQuantity>())
      {
         if (auto method = entity.get("MethodOfMeasurement"); !method.isNull())
            pset.set_attribute_value("MethodOfMeasurement", (std::string)method);
         pset.set_attribute_value("Quantities", std::vector<express::base>{});
      }
      else
      {
         pset.set_attribute_value("HasProperties", std::vector<express::base>{});
      }
      auto rel = m_Src.create<IfcSchema::IfcRelDefinesByProperties>();
      rel.setGlobalId(ifcopenshell::global_id());
      rel.setRelatedObjects(std::vector<IfcSchema::IfcObjectDefinition>{ object });
      rel.setRelatingPropertyDefinition(IfcSchema::IfcPropertySetDefinitionSelect(pset));
      return pset.as<IfcSchema::IfcPropertySetDefinition>();
   }

   std::pair<IndexType, IndexType> Merger::UpdateProperties(const std::string& key, PropertyOwner attach, const express::base& exported_pset, express::entity source_pset, const std::string& where, bool bShared)
   {
      IndexType nChanged = 0, nAdded = 0;
      std::string pset_name = exported_pset.as<express::entity>().get("Name");

      // the property list attribute of each kind of property container
      auto list_attribute = [](const express::base& container) -> std::string
         {
            if (container.as<IfcSchema::IfcElementQuantity>()) return "Quantities";
            if (container.as<IfcSchema::IfcMaterialProperties>()) return "Properties";
            return "HasProperties";
         };
      std::vector<express::base> exported_properties = exported_pset.as<express::entity>().get(list_attribute(exported_pset));
      std::vector<express::base> source_properties = source_pset.get(list_attribute(source_pset));
      bool bListChanged = false;

      for (const auto& exported_property : exported_properties)
      {
         std::string name = exported_property.as<express::entity>().get("Name");
         const auto* declaration = FindDeclaration(key, attach, pset_name, name);
         bool bTarget = declaration && declaration->target && !declaration->value;

         auto found = std::find_if(source_properties.begin(), source_properties.end(), [&](const express::base& p) {return (std::string)p.as<express::entity>().get("Name") == name; });
         if (found == source_properties.end())
         {
            source_properties.push_back(exported_property.as<IfcSchema::IfcPhysicalQuantity>() ? CopyQuantity(exported_property) : CopyProperty(exported_property));
            bListChanged = true;
            nAdded++;
            m_Result.log.push_back(where + ": added " + pset_name + "." + name);
            continue;
         }
         if (!bTarget)
            continue; // a constant or a placeholder: the source knows better

         // PGSuper's value, in the source property's unit
         express::base source_property = *found;
         express::base new_property;
         std::string old_text, new_text;
         if (auto sq = source_property.as<IfcSchema::IfcPhysicalSimpleQuantity>())
         {
            auto value_attribute = quantity_value_of(source_property.declaration().name());
            auto eq = exported_property.as<IfcSchema::IfcPhysicalSimpleQuantity>();
            if (!value_attribute || !eq || exported_property.declaration().name() != source_property.declaration().name())
               continue;
            double v = exported_property.as<express::entity>().get(value_attribute->first);
            double from = m_ExpUnits.GetConversionFactor(IfcSchema::IfcUnit(eq.Unit()), value_attribute->second);
            double to = m_SrcUnits.GetConversionFactor(IfcSchema::IfcUnit(sq.Unit()), value_attribute->second);
            double old_value = source_property.as<express::entity>().get(value_attribute->first);
            double new_value = v * from / to;
            if (IsEqual(old_value, new_value, 1.0e-9 * std::max(1.0, fabs(new_value))))
               continue;
            if (bShared || 1 < m_Src.get_total_inverses(source_property.id()))
            {
               new_property = m_Src.add_entity(exported_property);
               new_property.as<IfcSchema::IfcPhysicalSimpleQuantity>().setUnit(sq.Unit());
            }
            else
            {
               new_property = source_property;
            }
            new_property.set_attribute_value(value_attribute->first, new_value);
            old_text = std::to_string(old_value);
            new_text = std::to_string(new_value);
         }
         else
         {
            auto ss = source_property.as<IfcSchema::IfcPropertySingleValue>();
            auto es = exported_property.as<IfcSchema::IfcPropertySingleValue>();
            if (ss && es)
            {
               auto source_value = ss.NominalValue();
               auto exported_value = es.NominalValue();
               if (!exported_value)
                  continue;

               express::base value;
               bool bSameKind = source_value && unit_type_of_measure(source_value.declaration().name()) && unit_type_of_measure(source_value.declaration().name()) == unit_type_of_measure(exported_value.declaration().name());
               if (bSameKind)
               {
                  // a measure of the same kind: the source's type and unit
                  double from = m_ExpUnits.GetConversionFactor(es.Unit(), *unit_type_of_measure(exported_value.declaration().name()));
                  double to = m_SrcUnits.GetConversionFactor(ss.Unit(), *unit_type_of_measure(source_value.declaration().name()));
                  double v = *number_of(exported_value) * from / to;
                  auto old_v = number_of(source_value);
                  if (old_v && IsEqual(*old_v, v, 1.0e-9 * std::max(1.0, fabs(v))))
                     continue;
                  value = m_Src.create(&source_value.declaration());
                  value.set_attribute_value(0, v);
               }
               else if (source_value && source_value.declaration().name() == exported_value.declaration().name())
               {
                  if (value_text(source_value) == value_text(exported_value))
                     continue;
                  value = m_Src.add_entity(exported_value);
               }

               if (value)
               {
                  old_text = value_text(source_value);
                  new_text = value_text(value);
                  if (bShared || 1 < m_Src.get_total_inverses(source_property.id()))
                  {
                     auto copy = m_Src.create<IfcSchema::IfcPropertySingleValue>();
                     copy.setName(ss.Name());
                     copy.setSpecification(ss.Specification());
                     copy.setUnit(ss.Unit());
                     new_property = copy;
                  }
                  else
                  {
                     new_property = source_property;
                  }
                  new_property.as<IfcSchema::IfcPropertySingleValue>().setNominalValue(IfcSchema::IfcValue(value));
               }
               else
               {
                  // another kind of value (e.g. text where PGSuper writes a measure): PGSuper's property replaces it
                  new_property = CopyProperty(exported_property);
                  old_text = value_text(source_value);
                  new_text = value_text(new_property.as<IfcSchema::IfcPropertySingleValue>().NominalValue());
               }
            }
            else if (exported_property.declaration().name() == source_property.declaration().name() && exported_property.as<IfcSchema::IfcPropertyEnumeratedValue>())
            {
               auto ee = exported_property.as<IfcSchema::IfcPropertyEnumeratedValue>();
               auto se = source_property.as<IfcSchema::IfcPropertyEnumeratedValue>();
               auto text = [](const std::optional<std::vector<IfcSchema::IfcValue>>& values)
                  {
                     std::string t;
                     for (const auto& v : values.value_or(std::vector<IfcSchema::IfcValue>{}))
                        t += value_text(v) + ";";
                     return t;
                  };
               old_text = text(se.EnumerationValues());
               new_text = text(ee.EnumerationValues());
               if (old_text == new_text)
                  continue;
               std::vector<IfcSchema::IfcValue> values;
               for (const auto& v : ee.EnumerationValues().value_or(std::vector<IfcSchema::IfcValue>{}))
                  values.push_back(IfcSchema::IfcValue(m_Src.add_entity(v)));
               if (bShared || 1 < m_Src.get_total_inverses(source_property.id()))
               {
                  auto copy = m_Src.create<IfcSchema::IfcPropertyEnumeratedValue>();
                  copy.setName(se.Name());
                  copy.setSpecification(se.Specification());
                  copy.setEnumerationReference(se.EnumerationReference());
                  new_property = copy;
               }
               else
               {
                  new_property = source_property;
               }
               new_property.as<IfcSchema::IfcPropertyEnumeratedValue>().setEnumerationValues(values);
            }
            else
            {
               if (exported_property.declaration().name() == source_property.declaration().name())
                  continue; // e.g. list values: the source's are kept
               new_property = CopyProperty(exported_property);
               old_text = source_property.declaration().name();
               new_text = exported_property.declaration().name();
            }
         }

         if (new_property != source_property)
         {
            *found = new_property;
            bListChanged = true;
         }
         nChanged++;
         m_Result.log.push_back(where + ": " + pset_name + "." + name + " " + old_text + " -> " + new_text);
      }

      if (bListChanged)
         source_pset.set_attribute_value(list_attribute(source_pset), source_properties);
      return { nChanged, nAdded };
   }

   void Merger::UpdateElement(const Match& match)
   {
      auto source = match.source;
      std::string where = Describe(source);
      IndexType nChanged = 0, nAdded = 0;
      auto add = [&](std::pair<IndexType, IndexType> counts) { nChanged += counts.first; nAdded += counts.second; };

      // occurrence property sets and quantity sets
      for (const auto& rel : defined_by(match.exported))
      {
         auto exported_pset = rel.RelatingPropertyDefinition().as<IfcSchema::IfcPropertySetDefinition>();
         if (!exported_pset || !(exported_pset.as<IfcSchema::IfcPropertySet>() || exported_pset.as<IfcSchema::IfcElementQuantity>()))
            continue;
         std::string name = exported_pset.Name().value_or("");

         IfcSchema::IfcRelDefinesByProperties source_rel;
         auto source_pset = FindOccurrencePset(source, name, &source_rel);
         if (source_pset && source_pset.declaration().name() != exported_pset.declaration().name())
            continue; // e.g. a quantity set with a property set's name
         if (!source_pset)
            source_pset = NewOccurrencePset(source, exported_pset);
         else
            source_pset = OwnOccurrencePset(source, source_pset, source_rel);
         add(UpdateProperties(match.key, PropertyOwner::Occurrence, exported_pset, source_pset.as<express::entity>(), where, false));
      }

      if (auto exported = match.exported.as<IfcSchema::IfcObject>())
      {
         // type property sets: on the source's type, unless the source has no type or its type is shared by elements whose
         // values differ (then on the element, which overrides its type)
         for (const auto& type_rel : exported.IsTypedBy())
         {
            auto exported_type = type_rel.RelatingType();
            for (const auto& definition : exported_type.HasPropertySets().value_or(std::vector<IfcSchema::IfcPropertySetDefinition>{}))
            {
               auto exported_pset = definition.as<IfcSchema::IfcPropertySet>();
               if (!exported_pset)
                  continue;
               std::string name = exported_pset.Name().value_or("");

               IfcSchema::IfcTypeObject source_type;
               if (auto so = source.as<IfcSchema::IfcObject>(); so && !so.IsTypedBy().empty())
                  source_type = so.IsTypedBy().front().RelatingType();

               IfcSchema::IfcPropertySet type_pset;
               if (source_type)
               {
                  for (const auto& d : source_type.HasPropertySets().value_or(std::vector<IfcSchema::IfcPropertySetDefinition>{}))
                  {
                     if (d.Name().value_or("") == name && d.as<IfcSchema::IfcPropertySet>())
                        type_pset = d.as<IfcSchema::IfcPropertySet>();
                  }
               }

               // the values this element would give the type's property set
               std::string values;
               for (const auto& p : exported_pset.HasProperties())
               {
                  if (auto s = p.as<IfcSchema::IfcPropertySingleValue>())
                     values += p.Name() + "=" + value_text(s.NominalValue()) + ";";
               }
               bool bTypeIsShared = source_type && 1 < source_type.Types().size();
               auto shared_key = std::make_pair(source_type ? source_type.id() : 0u, name);
               auto assigned = m_SharedValues.find(shared_key);
               bool bConflict = source_type && assigned != m_SharedValues.end() && assigned->second != values;

               if (source_type && !bConflict)
               {
                  m_SharedValues[shared_key] = values;
                  if (!type_pset)
                  {
                     type_pset = m_Src.create<IfcSchema::IfcPropertySet>();
                     type_pset.setGlobalId(ifcopenshell::global_id());
                     type_pset.setName(name);
                     type_pset.setDescription(exported_pset.Description());
                     type_pset.setHasProperties(std::vector<IfcSchema::IfcProperty>{});
                     auto has = source_type.HasPropertySets().value_or(std::vector<IfcSchema::IfcPropertySetDefinition>{});
                     has.push_back(type_pset);
                     source_type.setHasPropertySets(has);
                  }
                  add(UpdateProperties(match.key, PropertyOwner::Type, exported_pset, type_pset.as<express::entity>(), where + " (its type " + Describe(source_type) + ")", bTypeIsShared));
               }
               else
               {
                  IfcSchema::IfcRelDefinesByProperties source_rel;
                  auto source_pset = FindOccurrencePset(source, name, &source_rel);
                  if (!source_pset)
                     source_pset = NewOccurrencePset(source, exported_pset);
                  else
                     source_pset = OwnOccurrencePset(source, source_pset, source_rel);
                  add(UpdateProperties(match.key, PropertyOwner::Type, exported_pset, source_pset.as<express::entity>(), where, false));
               }
            }
         }

         // material properties: on the source's material, which is copied for the element if it's shared by elements whose
         // values differ. An element without a material gets PGSuper's
         for (const auto& association : exported.HasAssociations())
         {
            auto material_rel = association.as<IfcSchema::IfcRelAssociatesMaterial>();
            if (!material_rel)
               continue;
            auto exported_material = material_rel.RelatingMaterial().as<IfcSchema::IfcMaterial>();
            if (!exported_material)
               continue;

            IfcSchema::IfcRelAssociatesMaterial source_material_rel;
            for (const auto& a : source.HasAssociations())
            {
               if (auto m = a.as<IfcSchema::IfcRelAssociatesMaterial>())
                  source_material_rel = m;
            }

            if (!source_material_rel)
            {
               auto material = m_Src.add_entity(exported_material).as<IfcSchema::IfcMaterial>();
               for (const auto& properties : exported_material.HasProperties())
               {
                  if (properties.Material() == exported_material)
                     m_Src.add_entity(properties);
               }
               auto rel = m_Src.create<IfcSchema::IfcRelAssociatesMaterial>();
               rel.setGlobalId(ifcopenshell::global_id());
               rel.setRelatedObjects(std::vector<IfcSchema::IfcDefinitionSelect>{ source });
               rel.setRelatingMaterial(material);
               m_Result.log.push_back(where + ": added PGSuper's material \"" + exported_material.Name() + "\"");
               nAdded++;
               continue;
            }

            auto source_material = source_material_rel.RelatingMaterial().as<IfcSchema::IfcMaterial>();
            if (!source_material)
            {
               m_Result.log.push_back(where + ": the material isn't a single IfcMaterial (" + source_material_rel.RelatingMaterial().declaration().name() + "), so its properties aren't changed");
               continue;
            }

            // the values this element would give the material
            std::string values;
            for (const auto& properties : exported_material.HasProperties())
            {
               for (const auto& p : properties.Properties())
               {
                  if (auto s = p.as<IfcSchema::IfcPropertySingleValue>())
                     values += properties.Name().value_or("") + "." + p.Name() + "=" + value_text(s.NominalValue()) + ";";
               }
            }
            auto shared_key = std::make_pair(source_material.id(), std::string("material"));
            auto assigned = m_SharedValues.find(shared_key);
            if (assigned != m_SharedValues.end() && assigned->second != values)
            {
               // another element gave the material other values: the element gets a copy of it, shared with the elements
               // that give it the same values
               auto& copy = m_MaterialCopies[std::make_pair(source_material.id(), values)];
               bool bNewCopy = !copy;
               if (bNewCopy)
               {
                  auto material = m_Src.create<IfcSchema::IfcMaterial>();
                  material.setName(source_material.Name());
                  material.setDescription(source_material.Description());
                  material.setCategory(source_material.Category());
                  for (const auto& properties : source_material.HasProperties())
                  {
                     auto pc = m_Src.create<IfcSchema::IfcMaterialProperties>();
                     pc.setName(properties.Name());
                     pc.setDescription(properties.Description());
                     pc.setProperties(properties.Properties()); // shared: a property that changes is replaced in the copy only
                     pc.setMaterial(material);
                  }
                  auto rel = m_Src.create<IfcSchema::IfcRelAssociatesMaterial>();
                  rel.setGlobalId(ifcopenshell::global_id());
                  rel.setRelatedObjects(std::vector<IfcSchema::IfcDefinitionSelect>{});
                  rel.setRelatingMaterial(material);
                  copy = material;
               }
               auto related = source_material_rel.RelatedObjects();
               related.erase(std::remove(related.begin(), related.end(), IfcSchema::IfcDefinitionSelect(source)), related.end());
               source_material_rel.setRelatedObjects(related);
               auto copy_rel = copy.AssociatedTo().front();
               auto copy_related = copy_rel.RelatedObjects();
               copy_related.push_back(source);
               copy_rel.setRelatedObjects(copy_related);
               m_Result.log.push_back(where + ": the material \"" + source_material.Name() + "\" is shared with elements whose PGSuper values differ, so the element " + (bNewCopy ? "gets a copy of it" : "shares the copy of it with the same values"));
               source_material = copy.as<IfcSchema::IfcMaterial>();
               shared_key = std::make_pair(copy.id(), std::string("material"));
            }
            m_SharedValues[shared_key] = values;

            for (const auto& exported_properties : exported_material.HasProperties())
            {
               std::string name = exported_properties.Name().value_or("");
               IfcSchema::IfcMaterialProperties source_properties;
               for (const auto& p : source_material.HasProperties())
               {
                  if (p.Name().value_or("") == name)
                     source_properties = p;
               }
               if (!source_properties)
               {
                  source_properties = m_Src.create<IfcSchema::IfcMaterialProperties>();
                  source_properties.setName(name);
                  source_properties.setDescription(exported_properties.Description());
                  source_properties.setProperties(std::vector<IfcSchema::IfcProperty>{});
                  source_properties.setMaterial(source_material);
               }
               add(UpdateProperties(match.key, PropertyOwner::Material, exported_properties, source_properties.as<express::entity>(), where + " (its material \"" + source_material.Name() + "\")", false));
            }
         }
      }

      UpdateClassifications(match);

      m_Result.nPropertiesChanged += nChanged;
      m_Result.nPropertiesAdded += nAdded;
      if (0 < nChanged + nAdded)
         m_Result.nUpdated++;
   }

   void Merger::UpdateClassifications(const Match& match)
   {
      auto reference_key = [](IfcSchema::IfcClassificationReference reference)
         {
            std::string system;
            if (auto c = reference.ReferencedSource().as<IfcSchema::IfcClassification>())
               system = c.Name();
            return system + "|" + reference.Identification().value_or("");
         };

      std::set<std::string> existing;
      for (const auto& association : match.source.HasAssociations())
      {
         if (auto rel = association.as<IfcSchema::IfcRelAssociatesClassification>())
         {
            if (auto reference = rel.RelatingClassification().as<IfcSchema::IfcClassificationReference>())
               existing.insert(reference_key(reference));
         }
      }

      for (const auto& association : match.exported.HasAssociations())
      {
         auto rel = association.as<IfcSchema::IfcRelAssociatesClassification>();
         if (!rel)
            continue;
         auto reference = rel.RelatingClassification().as<IfcSchema::IfcClassificationReference>();
         if (!reference || existing.count(reference_key(reference)))
            continue;

         // the source's reference with the same system and identification, else a copy of PGSuper's
         IfcSchema::IfcClassificationReference source_reference;
         for (const auto& r : m_Src.instances_by_type<IfcSchema::IfcClassificationReference>())
         {
            if (reference_key(r) == reference_key(reference))
            {
               source_reference = r;
               break;
            }
         }
         if (!source_reference)
         {
            source_reference = m_Src.add_entity(reference).as<IfcSchema::IfcClassificationReference>();
            MergeCopiedContextsAndClassifications(); // the copied system is the source's if it has one
         }

         IfcSchema::IfcRelAssociatesClassification source_rel;
         auto rels = source_reference.ClassificationRefForObjects();
         if (!rels.empty())
         {
            source_rel = rels.front();
            auto related = source_rel.RelatedObjects();
            related.push_back(match.source);
            source_rel.setRelatedObjects(related);
         }
         else
         {
            source_rel = m_Src.create<IfcSchema::IfcRelAssociatesClassification>();
            source_rel.setGlobalId(ifcopenshell::global_id());
            source_rel.setRelatedObjects(std::vector<IfcSchema::IfcDefinitionSelect>{ match.source });
            source_rel.setRelatingClassification(source_reference);
         }
         existing.insert(reference_key(reference));
         m_Result.nPropertiesAdded++;
         m_Result.log.push_back(Describe(match.source) + ": added the classification " + reference_key(reference));
      }
   }

   void Merger::Run()
   {
      std::ostringstream os;
      os << "Source model length unit: " << m_SrcMeters << " m; PGSuper's model: " << m_ExpMeters << " m.";
      Log(os.str());

      MatchElements();
      AddElements();
      for (const auto& match : m_Matches)
         UpdateElement(match);
   }
}

RoundTripMergeResult MergeIntoSourceModel(ifcopenshell::file& exported, const IfcExportIdentities& identities, const CRoundTripRegistry& registry, const CIfcMappingTable& table,
   const std::filesystem::path& source_path, const std::filesystem::path& result_path)
{
   USES_CONVERSION;
   RoundTripMergeResult result;

   ifcopenshell::file source(source_path.string());
   if (!source.good())
      throw std::runtime_error("The IFC model the project was imported from can't be read: " + source_path.string());
   if (source.schema() != exported.schema())
      throw std::runtime_error("The IFC model the project was imported from is " + source.schema()->name() + ", but the export is " + exported.schema()->name() + ".");

   Merger merger(exported, source, identities, registry, table, result);
   merger.Run();

   std::ofstream ofs(result_path);
   if (!ofs)
      throw std::runtime_error("The IFC file can't be written: " + result_path.string());
   ofs << source;
   return result;
}
