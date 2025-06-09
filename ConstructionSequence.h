///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright © 1999-2025  Washington State Department of Transportation
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

#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>

// Per AASHTO LRFD 2.5.3, "When the designer has assumed a particular sequence of construction in order to induce
// certain stresses under dead load, that sequence shall be defined in the contract documents."
// 
// Also, per 5.9.4.5 "Detensioning of temporary strands shall be shown in the construction sequence and typically
// occurs after the girders are securely braced and before construction of intermediate concrete diaphragm, if applicable."
//
// The CreateWorkPlan function defines the assumed construction sequence as an IfcWorkPlan
// 
// See discussion at https://community.osarch.org/discussion/comment/23500
//
// IfcProject - IfcRelDeclares - IfcWorkPlan - IfcRelAggregates - IfcWorkSchedule - IfcRelAssignsToControl - IfcTask (summary) - IfcRelNests - [Subtasks IfcTask, IfcTask, etc]
//                                                              |                                                 |
//                                                              |                                        IfcRelAssignsToProduct
//                                                              |                                                 |
//                                                              |                                             IfcBridge (bridge #1)
//                                                              |
//                                                              |
//                                                              - IfcWorkSchedule - IfcRelAssignsToControl - IfcTask (summary) - IfcRelNests - [Subtasks IfcTask, IfcTask, etc]
//                                                                                                                |
//                                                                                                       IfcRelAssignsToProduct
//                                                                                                                |
//                                                                                                            IfcBridge (bridge #2)

// This name is the lookup key for the assumed construction sequence work schedule.
// Using a global variable so it can be used in the definition and the look up.
static std::string gs_work_schedule_name = std::string("Assumed Construction Sequence");

// this function generated from Microsoft Copilot
std::string getCurrentISO8601Time() {
   // Get the current time
   auto now = std::chrono::system_clock::now();
   std::time_t now_time = std::chrono::system_clock::to_time_t(now);

   // Convert to tm structure
   std::tm tm = *std::gmtime(&now_time);

   // Format the time to ISO 8601
   std::ostringstream oss;
   oss << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
   return oss.str();
}



template <typename Schema>
void CreateAssumedConstructionSequence(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcModelBuilderOptions& options)
{
   // only doing the PGSuper assumed construction sequence for now. For time-step method (PGSplice) it is much more complex
   USES_CONVERSION;

   GET_IFACE2_NOCHECK(pBroker, ILossParameters, pLossParams);

   if (!options.include_work_plan || pLossParams->GetLossMethod() == PrestressLossCriteria::LossMethodType::TIME_STEP)
      return;


   auto now = getCurrentISO8601Time();

   // Define the work plan

   auto work_plan = new Schema::IfcWorkPlan(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Bridge Project Work Plan"), // Name
      std::string("Work schedules for construction"), // Description
      boost::none, // ObjectType
      boost::none, // Identification
      now, // CreationData
      boost::none, // Creators
      std::string("Defines assumed plans and schedules required per AASHTO specifications"), // Purpose
      boost::none, // Duration
      boost::none, // TotalFloat
      now, // StartTime (we don't know the start time, so just use the creation time)
      boost::none, // FinishTime
      Schema::IfcWorkPlanTypeEnum::IfcWorkPlanType_PLANNED
   );

   // Declare the work plan in the project
   auto project = file.getSingle<typename Schema::IfcProject>();
   file.addRelatedObject<typename Schema::IfcRelDeclares>(project, work_plan);
   
   // Define the assumed construction sequence as a work schedule
   auto work_schedule = new Schema::IfcWorkSchedule(
      IfcParse::IfcGlobalId(),
      nullptr,
      gs_work_schedule_name, // Name
      boost::none, // Description
      boost::none, // ObjectType
      boost::none, // Identification
      now, // CreationDate
      boost::none, // Creators
      std::string("Satisfies the requirements of AASHTO LRFD BDS 2.5.3"), // Purpose
      boost::none, // Duration
      boost::none, // TotalFloat
      now, // StartTime
      boost::none, // FinishTime
      Schema::IfcWorkScheduleTypeEnum::IfcWorkScheduleType_PLANNED // PredefinedType
   );


   // Aggregate work schedules with the work plan
   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr work_schedules(new aggregate_of<typename Schema::IfcObjectDefinition>());
   work_schedules->push(work_schedule);

   auto rel_aggregates = new Schema::IfcRelAggregates(
      IfcParse::IfcGlobalId(),
      nullptr,
      boost::none, // Name
      boost::none, // Description, 
      work_plan,
      work_schedules
   );

   file.addEntity(rel_aggregates);

   // Define the tasks associated with the work schedule

   // Define summary task for bridge construction
   // See https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/lexical/IfcTask.htm Object Nesting concept for
   // discussion about the Identification attribute
   auto summary_task = new Schema::IfcTask(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Construct Bridge"), // Name
      std::string("Assumed construction sequence"), // Description
      boost::none, // ObjectType
      std::string("0"), // Identification
      std::string("Provided per requirements of AASHTO LRFD BDS 2.5.3"), // Purpose
      boost::none, // Status
      boost::none, // WorkMethod
      true, // IsMilestone
      boost::none, // Priority, 
      nullptr, // TaskTime, 
      Schema::IfcTaskTypeEnum::IfcTaskType_CONSTRUCTION
   );

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr tasks(new aggregate_of<typename Schema::IfcObjectDefinition>());
   tasks->push(summary_task);

   // Assign the tasks of the construction sequence to the work schedule
   auto rel_assigns_to_control = new Schema::IfcRelAssignsToControl(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Assumed Construction Sequence Work Schedule Tasks"), // Name
      boost::none, // Description, 
      tasks, // RelatedObjects
      boost::none, // RelatedObjectsType
      work_schedule // RelatingControl
   );

   file.addEntity(rel_assigns_to_control);

   // Assign the task of construction to the bridge (bridge is the output of the task)
   auto bridge = file.getSingle<typename Schema::IfcBridge>();
   auto rel_assigns_to_product = new Schema::IfcRelAssignsToProduct(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Bridge is the output product of the construction task"), // Name
      boost::none, // Description
      tasks, // RelatedObjects
      boost::none, // RelatedObjectsType
      bridge // RelatingProduct
   );
   file.addEntity(rel_assigns_to_product);

   //
   // This is just a concept, but if the construction sequence is document based information, the document
   // can be associated with the summary task
   //

   // define document reference
   GET_IFACE2(pBroker, IBridge, pBridge);
   auto nSpans = pBridge->GetSpanCount();
   auto document_reference = new Schema::IfcDocumentReference(
      nSpans == 1 ? std::string("https://wsdot.wa.gov/publications/fulltext/Bridge/Web_BSD/5.6_A2_1.PDF") 
                  : std::string("https://wsdot.wa.gov/publications/fulltext/Bridge/Web_BSD/5.6_A2_2.PDF"),
      nSpans == 1 ? std::string("5.6 A2 1") : std::string("5.6 A2 2"), // Identification
      std::string("Assumed construction sequence"), // Name
      boost::none, // Description 
      nullptr);

   typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr tasks2(new aggregate_of<typename Schema::IfcDefinitionSelect>());
   tasks2->push(summary_task);

   // associate document with task
   auto rel_associates_document = new Schema::IfcRelAssociatesDocument(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, tasks2, document_reference);
   file.addEntity(rel_associates_document);

   //
   // Define sub-tasks of the summary task
   //
   
   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr sub_tasks(new aggregate_of<typename Schema::IfcObjectDefinition>());

   auto sub_task1 = CreateStage1Tasks(file, pBroker, options);
   auto sub_task2 = CreateStage2Tasks(file, pBroker, options);
   auto sub_task3 = CreateStage3Tasks(file, pBroker, options);
   auto sub_task4 = CreateStage4Tasks(file, pBroker, options);

   sub_tasks->push(sub_task1);
   sub_tasks->push(sub_task2);
   sub_tasks->push(sub_task3);
   sub_tasks->push(sub_task4);

   // nest the subtasks into the summary task
   auto rel_nests = new Schema::IfcRelNests(
      IfcParse::IfcGlobalId(),
      nullptr,
      boost::none, // Name
      boost::none, // Description
      summary_task, // RelatingObject
      sub_tasks // RelatedObjects
   );
   file.addEntity(rel_nests);

   // Assign relative sequence to subtasks
   auto iter = std::next(sub_tasks->begin());
   auto end = sub_tasks->end();
   for (; iter != end; iter++)
   {
      auto object_definition = *(std::prev(iter));
      auto relating_process = object_definition->as<typename Schema::IfcProcess>();
      object_definition = *iter;
      auto related_process = object_definition->as<typename Schema::IfcProcess>();

      auto rel_sequence = new Ifc4x3_add2::IfcRelSequence(
            IfcParse::IfcGlobalId(),
            nullptr,
            boost::none, // Name
            boost::none, // Description
            relating_process,
            related_process,
            nullptr, // TimeLag
            Schema::IfcSequenceEnum::IfcSequence_FINISH_START, // SequenceType
            boost::none // UserDefinedSequenceType
         );
      file.addEntity(rel_sequence);
   }

}

template <typename Schema>
typename Schema::IfcTask* CreateStage1Tasks(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcModelBuilderOptions& options)
{
   auto task = new Schema::IfcTask(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Stage 1 - Set girders in place"), // Name
      boost::none, // Description
      boost::none, // ObjectType
      std::string("1"), // Identification
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

   auto construct_piers = new Schema::IfcTask(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Construct piers and abutments"), // Name
      boost::none, // Description
      boost::none, // ObjectType
      std::string("1.1"), // Identification
      boost::none, // LongDescription
      boost::none, // Status
      boost::none, // WorkMethod
      true, // IsMilestone
      boost::none, // Priority, 
      nullptr, // TaskTime, 
      Schema::IfcTaskTypeEnum::IfcTaskType_CONSTRUCTION
   );

   auto install_bearings = new Schema::IfcTask(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Install elastomeric bearing pads"), // Name
      boost::none, // Description
      boost::none, // ObjectType
      std::string("1.2"), // Identification
      boost::none, // LongDescription
      boost::none, // Status
      boost::none, // WorkMethod
      true, // IsMilestone
      boost::none, // Priority, 
      nullptr, // TaskTime, 
      Schema::IfcTaskTypeEnum::IfcTaskType_INSTALLATION
   );

   auto set_girders = new Schema::IfcTask(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Set and brace girders"), // Name
      boost::none, // Description
      boost::none, // ObjectType
      std::string("1.3"), // Identification
      std::string("Install temporary bracing for erection in accordance with Std. Spec. Section 6-02.3(17)F4."), // LongDescription
      boost::none, // Status
      boost::none, // WorkMethod
      true, // IsMilestone
      boost::none, // Priority, 
      nullptr, // TaskTime, 
      Schema::IfcTaskTypeEnum::IfcTaskType_INSTALLATION
   );

   subtasks->push(construct_piers);
   subtasks->push(install_bearings);
   subtasks->push(set_girders);

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
         std::string("1.4"), // Identification
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
         std::string("1.4.1"), // Identification
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
         std::string("1.4.2"), // Identification
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
         std::string("1.4.3"), // Identification
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

      // connnect sequence of tasks
      auto iter = std::next(remove_temp_strands_subtasks->begin());
      auto end = remove_temp_strands_subtasks->end();
      for (; iter != end; iter++)
      {
         auto object_definition = *(std::prev(iter));
         auto relating_process = object_definition->as<typename Schema::IfcProcess>();
         object_definition = *iter;
         auto related_process = object_definition->as<typename Schema::IfcProcess>();

         auto rel_sequence = new Ifc4x3_add2::IfcRelSequence(
            IfcParse::IfcGlobalId(),
            nullptr,
            boost::none, // Name
            boost::none, // Description
            relating_process,
            related_process,
            nullptr, // TimeLag
            Schema::IfcSequenceEnum::IfcSequence_FINISH_START, // SequenceType
            boost::none // UserDefinedSequenceType
         );
         file.addEntity(rel_sequence);
      }

   }

   auto rel_nests = new Schema::IfcRelNests(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, task, subtasks);
   file.addEntity(rel_nests);

   // connect sequence of tasks
   auto iter = std::next(subtasks->begin());
   auto end = subtasks->end();
   for (; iter != end; iter++)
   {
      auto object_definition = *(std::prev(iter));
      auto relating_process = object_definition->as<typename Schema::IfcProcess>();
      object_definition = *iter;
      auto related_process = object_definition->as<typename Schema::IfcProcess>();

      auto rel_sequence = new Ifc4x3_add2::IfcRelSequence(
         IfcParse::IfcGlobalId(),
         nullptr,
         boost::none, // Name
         boost::none, // Description
         relating_process,
         related_process,
         nullptr, // TimeLag
         Schema::IfcSequenceEnum::IfcSequence_FINISH_START, // SequenceType
         boost::none // UserDefinedSequenceType
      );
      file.addEntity(rel_sequence);
   }

   return task;
}

template <typename Schema>
typename Schema::IfcTask* CreateStage2Tasks(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcModelBuilderOptions& options)
{
   auto task = new Schema::IfcTask(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Stage 2 - Cast diaphragms and place bridge deck reinforcement"), // Name
      boost::none, // Description
      boost::none, // ObjectType
      std::string("2"), // Identification
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
      std::string("Install temporary bracing"), // Name
      boost::none, // Description
      boost::none, // ObjectType
      std::string("2.1"), // Identification
      std::string("Install temporary bracing for diaphragm and deck concrete placement in accordance with Std. Spec. Section 6-02.3(17)F5."), // LongDescription
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
      std::string("Construct deck"), // Name
      boost::none, // Description
      boost::none, // ObjectType
      std::string("2.2"), // Identification
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


   auto iter = std::next(subtasks->begin());
   auto end = subtasks->end();
   for (; iter != end; iter++)
   {
      auto object_definition = *(std::prev(iter));
      auto relating_process = object_definition->as<typename Schema::IfcProcess>();
      object_definition = *iter;
      auto related_process = object_definition->as<typename Schema::IfcProcess>();

      auto rel_sequence = new Ifc4x3_add2::IfcRelSequence(
         IfcParse::IfcGlobalId(),
         nullptr,
         boost::none, // Name
         boost::none, // Description
         relating_process,
         related_process,
         nullptr, // TimeLag
         Schema::IfcSequenceEnum::IfcSequence_FINISH_START, // SequenceType
         boost::none // UserDefinedSequenceType
      );
      file.addEntity(rel_sequence);
   }

   return task;
}

template <typename Schema>
typename Schema::IfcTask* CreateStage3Tasks(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcModelBuilderOptions& options)
{
   auto task = new Schema::IfcTask(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Stage 3 - Cast bridge deck"), // Name
      boost::none, // Description
      boost::none, // ObjectType
      std::string("3"), // Identification
      std::string("Cast bridge deck when diaphragm concrete compressive strength has reached 3000 psi (minimum)."), // LongDescription
      boost::none, // Status
      boost::none, // WorkMethod
      true, // IsMilestone
      boost::none, // Priority, 
      nullptr, // TaskTime, 
      Schema::IfcTaskTypeEnum::IfcTaskType_CONSTRUCTION
   );

   return task;
}


template <typename Schema>
typename Schema::IfcTask* CreateStage4Tasks(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcModelBuilderOptions& options)
{
   auto task = new Schema::IfcTask(
      IfcParse::IfcGlobalId(),
      nullptr,
      std::string("Stage 4 - Cast traffic barriers"), // Name
      std::string("Cast traffic barrier"), // Description
      boost::none, // ObjectType
      std::string("4"), // Identification
      std::string("Cast traffic barrier after the deck concrete compressive strength has reached 3000 psi (minimum)"), // LongDescription
      boost::none, // Status
      boost::none, // WorkMethod
      true, // IsMilestone
      boost::none, // Priority, 
      nullptr, // TaskTime, 
      Schema::IfcTaskTypeEnum::IfcTaskType_CONSTRUCTION
   );

   return task;
}


template <typename Schema>
typename Schema::IfcTask* GetConstructBridgeSummaryTask(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcModelBuilderOptions& options)
{
   // the summary IfcTask is assigned to the IfcWorkSchedule control
   auto assignments = file.instances_by_type<typename Schema::IfcRelAssignsToControl>();
   for (auto& assignment : *assignments)
   {
      auto relating_control = assignment->RelatingControl()->as<typename Schema::IfcWorkSchedule>();
      if (relating_control && relating_control->Name() == gs_work_schedule_name) // is this an IfcWorkSchedule and does it have the right name?
      {
         auto object_definitions = assignment->RelatedObjects(); // IfcTask are the related objects
         auto object_definition = *(object_definitions->begin()); // We know there is only the one, so get it
         auto task = object_definition->as<typename Schema::IfcTask>(); // Recover its type
         return task;
      }
   }
   return nullptr;
}

template <typename Schema>
typename Schema::IfcTask* GetStageTask(IndexType idx,IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcModelBuilderOptions& options)
{
   auto summary_task = GetConstructBridgeSummaryTask(file, pBroker, options);
   auto nested_by = summary_task->IsNestedBy();
   auto rel_nests = *(nested_by->begin());
   auto related_objects = rel_nests->RelatedObjects();
   auto iter = related_objects->begin();
   iter += idx;
   auto object_definition = *iter;
   auto task = object_definition->as<typename Schema::IfcTask>();
   return task;
}

template <typename Schema>
typename Schema::IfcTask* GetStage1Task(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcModelBuilderOptions& options)
{
   return GetStageTask(0, file, pBroker, options);
}

template <typename Schema>
typename Schema::IfcTask* GetStage2Task(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcModelBuilderOptions& options)
{
   return GetStageTask(1, file, pBroker, options);
}

template <typename Schema>
typename Schema::IfcTask* GetStage3Task(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcModelBuilderOptions& options)
{
   return GetStageTask(2, file, pBroker, options);
}

template <typename Schema>
typename Schema::IfcTask* GetStage4Task(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcModelBuilderOptions& options)
{
   return GetStageTask(3, file, pBroker, options);
}

template <typename Schema>
void AssignTaskToProduct(typename Schema::IfcTask* task,typename Schema::IfcProduct* product,IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcModelBuilderOptions& options)
{
   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr tasks(new aggregate_of<typename Schema::IfcObjectDefinition>());
   tasks->push(task);
   auto rel_assigns_to_product = new Schema::IfcRelAssignsToProduct(
      IfcParse::IfcGlobalId(),
      nullptr,
      boost::none, // Name
      boost::none, // Description
      tasks, // RelatedObjects
      boost::none, // RelatedObjectsType
      product // RelatingProduct
   );
   file.addEntity(rel_assigns_to_product);
}