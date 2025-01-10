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

// Per AASHTO LRFD 2.5.3, "When the designer has assumed a particular sequence of constructon in order to induce
// certain stresses under dead load, that sequence shall be defined in the contract documents."
// 
// Also, per 5.9.4.5 "Detensioning of temporary strands shall be shown in the construction sequence and typically
// occurs after the girders are securely braced and before construction of intermediate concrete diaphragm, if applicable."
//
// The CreateWorkPlan function defines the assumed construction sequence as an IfcWorkPlan

template <typename Schema>
void CreateWorkPlan(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CIfcModelBuilderOptions& options)
{
   USES_CONVERSION;

   GET_IFACE2(pBroker, ILossParameters, pLossParams);
   // only doing the PGSuper assumed construction sequence for now. For time-step method (PGSplice) it is much more complex

   if (!options.include_work_plan || pLossParams->GetLossMethod() == PrestressLossCriteria::LossMethodType::TIME_STEP)
      return;

   auto work_plan = new Schema::IfcWorkPlan(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Assumed Construction Sequence"), // Name
      boost::none, // Description
      boost::none, // ObjectType
      boost::none, // Identification
      std::string("Unknown"), // CreationData
      boost::none, // Creators
      std::string("Satisfies the requirements of AASHTO LRFD BDS 2.5.3"), // Purpose
      boost::none, // Duration
      boost::none, // TotalFloat
      std::string("Unknown"), // StartTime
      boost::none, // FinishTime
      Schema::IfcWorkPlanTypeEnum::IfcWorkPlanType_PLANNED
   );

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr tasks(new aggregate_of<typename Schema::IfcObjectDefinition>());

   auto task1 = CreateStage1Tasks(file, pBroker, options);
   auto task2 = CreateStage2Tasks(file, pBroker, options);

   auto task3 = new Schema::IfcTask(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Stage 3"), // Name
      std::string("Cast bridge deck"), // Description
      boost::none, // ObjectType
      boost::none, // Identification
      std::string("Cast bridge deck when diaphragm concrete compressive strength has reached 3000 psi (minimum)."), // LongDescription
      boost::none, // Status
      boost::none, // WorkMethod
      true, // IsMilestone
      boost::none, // Priority, 
      nullptr, // TaskTime, 
      Schema::IfcTaskTypeEnum::IfcTaskType_CONSTRUCTION
   );

   auto task4 = new Schema::IfcTask(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Stage 4"), // Name
      std::string("Cast traffic barrier"), // Description
      boost::none, // ObjectType
      boost::none, // Identification
      std::string("Cast traffic barrier after the deck concrete compressive strength has reached 3000 psi (minimum)"), // LongDescription
      boost::none, // Status
      boost::none, // WorkMethod
      true, // IsMilestone
      boost::none, // Priority, 
      nullptr, // TaskTime, 
      Schema::IfcTaskTypeEnum::IfcTaskType_CONSTRUCTION
   );

   tasks->push(task1);
   tasks->push(task2);
   tasks->push(task3);
   tasks->push(task4);

   auto rel_assigns_to_control = new Ifc4x3_add2::IfcRelAssignsToControl(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Construction Sequence Tasks"), // Name
      boost::none, // Description, 
      tasks, // RelatedObjects
      boost::none, // RelatedObjectsType
      work_plan
   );

   file.addEntity(rel_assigns_to_control);



   auto project = file.getSingle<typename Schema::IfcProject>();
   auto rel_declares_instances = file.instances_by_type<typename Schema::IfcRelDeclares>();
   if (rel_declares_instances->size() == 0)
   {
      typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr related_definitions(new aggregate_of<typename Schema::IfcDefinitionSelect>());
      related_definitions->push(work_plan);

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
            related_definitions->push(work_plan);
         }
      }
   }
}

template <typename Schema>
typename Schema::IfcTask* CreateStage1Tasks(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CIfcModelBuilderOptions& options)
{
   auto task = new Schema::IfcTask(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Stage 1"), // Name
      std::string("Set girders in place"), // Description
      boost::none, // ObjectType
      boost::none, // Identification
      boost::none, // LongDescription
      boost::none, // Status
      boost::none, // WorkMethod
      true, // IsMilestone
      boost::none, // Priority, 
      nullptr, // TaskTime, 
      Schema::IfcTaskTypeEnum::IfcTaskType_MOVE
   );

   // https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/lexical/IfcTask.htm
   // Model the individual steps of girder erection step as sub-tasks of the parent task
   // Each sub-task is IfcRelNests with the parent task.
   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr subtasks(new aggregate_of<typename Schema::IfcObjectDefinition>());

   auto install_bracing_task = new Schema::IfcTask(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Brace girders"), // Name
      boost::none, // Description
      boost::none, // ObjectType
      boost::none, // Identification
      std::string("Install temporary bracing for erection in accordance with Std. Spec. Section 6-02.3(17)F4."), // LongDescription
      boost::none, // Status
      boost::none, // WorkMethod
      true, // IsMilestone
      boost::none, // Priority, 
      nullptr, // TaskTime, 
      Schema::IfcTaskTypeEnum::IfcTaskType_INSTALLATION
   );

   subtasks->push(install_bracing_task);

   bool bHasTempStrands = false;
   GET_IFACE2(pBroker, IStrandGeometry, pStrandGeom);
   GET_IFACE2(pBroker, IBridge, pBridge);
   GroupIndexType nGroups = pBridge->GetGirderGroupCount();
   for (GroupIndexType grpIdx = 0; grpIdx < nGroups; grpIdx++)
   {
      GirderIndexType nGirders = pBridge->GetGirderCount(grpIdx);
      for (GirderIndexType gdrIdx = 0; gdrIdx < nGirders; gdrIdx++)
      {
         SegmentIndexType nSegments = pBridge->GetSegmentCount(grpIdx, gdrIdx);
         for (SegmentIndexType segIdx = 0; segIdx < nSegments; segIdx++)
         {
            StrandIndexType nTempStrands = pStrandGeom->GetStrandCount(CSegmentKey(grpIdx, gdrIdx, segIdx), pgsTypes::Temporary);
            if (0 < nTempStrands)
            {
               bHasTempStrands = true;
               break;
            }
         }
         if (bHasTempStrands) break;
      }
      if (bHasTempStrands) break;
   }

   if (bHasTempStrands)
   {
      auto remove_temp_strands_task = new Schema::IfcTask(
         IfcParse::IfcGlobalId(),
         nullptr,
         std::string("Remove temporary strands"), // Name
         boost::none, // Description
         boost::none, // ObjectType
         boost::none, // Identification
         std::string("Satisfies the requirements of AASHTO LRFD BDS 5.9.4.5"), // LongDescription
         boost::none, // Status
         boost::none, // WorkMethod
         true, // IsMilestone
         boost::none, // Priority, 
         nullptr, // TaskTime, 
         Schema::IfcTaskTypeEnum::IfcTaskType_REMOVAL
      );

      subtasks->push(remove_temp_strands_task);

      typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr remove_temp_strands_subtasks(new aggregate_of<typename Schema::IfcObjectDefinition>());

      auto remove_poly = new Schema::IfcTask(
         IfcParse::IfcGlobalId(),
         nullptr,
         std::string("Remove polystyrene"), // Name
         boost::none, // Description
         boost::none, // ObjectType
         boost::none, // Identification
         std::string("Removed expanded polystyrene in blockouts in top flange of girders. Once the expanded polystyrene has been removed from the strand detensioning blockout, prevent moisture from entering the blockout until the temporary top strand is cut and the blockout filled with grout."), // LongDescription
         boost::none, // Status
         boost::none, // WorkMethod
         true, // IsMilestone
         boost::none, // Priority, 
         nullptr, // TaskTime, 
         Schema::IfcTaskTypeEnum::IfcTaskType_REMOVAL
      );

      auto cut_strands = new Schema::IfcTask(
         IfcParse::IfcGlobalId(),
         nullptr,
         std::string("Cut strands"), // Name
         boost::none, // Description
         boost::none, // ObjectType
         boost::none, // Identification
         std::string("Cut strands in blockouts. Strands may be cut by using a cutting torch and moving the flame back and forth over the length of the exposed strand to let individual wires break one at a time to lessen the shock to the girder. Strands shall be released in a symmetrical manner about the girder centerline starting with those furthest from the centerline and working inwards. For post-tensioned temporary top strands, actively restrain the strand chucks at the girder ends during cutting."), // LongDescription
         boost::none, // Status
         boost::none, // WorkMethod
         true, // IsMilestone
         boost::none, // Priority, 
         nullptr, // TaskTime, 
         Schema::IfcTaskTypeEnum::IfcTaskType_OPERATION
      );

      auto fill_blockouts = new Schema::IfcTask(
         IfcParse::IfcGlobalId(),
         nullptr,
         std::string("Fill blockouts"), // Name
         boost::none, // Description
         boost::none, // ObjectType
         boost::none, // Identification
         std::string("Within 24 hours of cutting the temporary strands, fill the blockouts with a grout conforming to Std. Spec. 9-20.3(2). Remove all moisture in blockouts prior to filling them with grout."), // LongDescription
         boost::none, // Status
         boost::none, // WorkMethod
         true, // IsMilestone
         boost::none, // Priority, 
         nullptr, // TaskTime, 
         Schema::IfcTaskTypeEnum::IfcTaskType_INSTALLATION
      );

      remove_temp_strands_subtasks->push(remove_poly);
      remove_temp_strands_subtasks->push(cut_strands);
      remove_temp_strands_subtasks->push(fill_blockouts);

      auto rel_nests = new Schema::IfcRelNests(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, remove_temp_strands_task, remove_temp_strands_subtasks);
      file.addEntity(rel_nests);
   }


   auto rel_nests = new Schema::IfcRelNests(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, task, subtasks);
   file.addEntity(rel_nests);

   return task;
}

template <typename Schema>
typename Schema::IfcTask* CreateStage2Tasks(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CIfcModelBuilderOptions& options)
{
   auto task = new Schema::IfcTask(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Stage 2"), // Name
      std::string("Cast diaphragms and place bridge deck reinforcement"), // Description
      boost::none, // ObjectType
      boost::none, // Identification
      boost::none, // LongDescription
      boost::none, // Status
      boost::none, // WorkMethod
      true, // IsMilestone
      boost::none, // Priority, 
      nullptr, // TaskTime, 
      Schema::IfcTaskTypeEnum::IfcTaskType_CONSTRUCTION
   );

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr subtasks(new aggregate_of<typename Schema::IfcObjectDefinition>());

   auto install_bracing = new Schema::IfcTask(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Diaphragm and deck placement temporary bracing"), // Name
      boost::none, // Description
      boost::none, // ObjectType
      boost::none, // Identification
      std::string("Install temporary bracing for diaphragm and deck placement in accordance with Std. Spec. Section 6-02.3(17)F5."), // LongDescription
      boost::none, // Status
      boost::none, // WorkMethod
      true, // IsMilestone
      boost::none, // Priority, 
      nullptr, // TaskTime, 
      Schema::IfcTaskTypeEnum::IfcTaskType_INSTALLATION
   );


   auto place_deck_reinforcement = new Schema::IfcTask(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Deck reinforcement"), // Name
      boost::none, // Description
      boost::none, // ObjectType
      boost::none, // Identification
      std::string("Form deck and place bridge deck reinforcement after casting diaphragms."), // LongDescription
      boost::none, // Status
      boost::none, // WorkMethod
      true, // IsMilestone
      boost::none, // Priority, 
      nullptr, // TaskTime, 
      Schema::IfcTaskTypeEnum::IfcTaskType_CONSTRUCTION
   );

   subtasks->push(install_bracing);
   subtasks->push(place_deck_reinforcement);



   auto rel_nests = new Schema::IfcRelNests(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, task, subtasks);
   file.addEntity(rel_nests);

   return task;
}