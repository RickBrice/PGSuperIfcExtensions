#pragma once

template <typename Schema>
void ProjectDeclares(IfcHierarchyHelper<Schema>& file, typename Schema::IfcDefinitionSelect* related_definition)
{
   // Declare work plan and work schedule in the project
   auto project = file.getSingle<typename Schema::IfcProject>();
   auto rel_declares_instances = file.instances_by_type<typename Schema::IfcRelDeclares>();
   if (rel_declares_instances->size() == 0)
   {
      typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr related_definitions(new aggregate_of<typename Schema::IfcDefinitionSelect>());
      related_definitions->push(related_definition);

      auto rel_declares = new Schema::IfcRelDeclares(
         IfcParse::IfcGlobalId(),
         nullptr,
         boost::none,
         boost::none,
         project,
         related_definitions);

      file.addEntity(rel_declares);
   }
   else
   {
      for (auto& rel_declares : *rel_declares_instances)
      {
         if (rel_declares->RelatingContext()->as<typename Schema::IfcProject>())
         {
            auto related_definitions = rel_declares->RelatedDefinitions();
            related_definitions->push(related_definition);
         }
      }
   }
}
