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

// Girder cross sections measured from the model, and compared with the sections of girder library entries
// (devdocs/GirderImportPlan.md, G1). Sections are in section coordinates: x horizontal, y up, origin at the
// top center of the section.

#include <Eigen/Dense>
#include <optional>
#include <string>
#include <vector>

struct Mesh;
class GirderLibraryEntry;
interface IShape;

// The outline of a section: its solid parts and their voids (e.g. the pieces of a composite library section). Solid
// parts don't overlap, and voids are inside them
struct SectionOutline
{
   using Loop = std::vector<Eigen::Vector2d>;
   std::vector<Loop> solids;
   std::vector<Loop> voids;
};

// The measures a section is compared by
struct SectionMeasures
{
   double height = 0;
   double area = 0;
   double centroid_y = 0; // from the top (negative)
   double ix = 0;         // about the horizontal axis through the centroid
   double width = 0;      // overall width
   std::vector<double> width_profile; // material width at each WIDTH_PROFILE_STEP below the top
};

constexpr double WIDTH_PROFILE_STEP = 0.003175; // 1/8 in

// The tolerances for a section to match another. The same for every girder family for now
// (devdocs/GirderImportPlan.md, decision 2), but kept as data so they can differ by family later
struct SectionTolerances
{
   double length = 0.003175; // 1/8 in: height and widths at each depth (and the depth where a width changes)
   double ratio = 0.01;      // area and moment of inertia
};

// The section of a girder mesh cut by a vertical plane normal to the girder axis (start to end, in plan) at a
// fraction of its length. Returns nullopt if the plane doesn't cut a closed outline
std::optional<SectionOutline> CutGirderSection(const Mesh& mesh, const Eigen::Vector2d& start, const Eigen::Vector2d& end, double fraction);

// The outline of a WBFL shape (e.g. a girder library entry's section), moved to section coordinates
SectionOutline GetShapeOutline(IShape* pShape);

// The outline of a girder library entry's section, or nullopt if its beam factory can't make one
std::optional<SectionOutline> GetLibraryEntryOutline(const GirderLibraryEntry* pEntry);

SectionMeasures MeasureSection(const SectionOutline& outline);

// The differences between two sections that are larger than the tolerances, as text (e.g. "height 1.250 in"). Empty
// if the sections match
std::vector<std::string> CompareSections(const SectionMeasures& model, const SectionMeasures& library, const SectionTolerances& tolerances);

// A number that grows with the difference between two sections, to rank the library entries
double SectionDifference(const SectionMeasures& model, const SectionMeasures& library);

// e.g. "height 72.000 in, area 913.4 in^2, I 545,894 in^4, top width 49.000 in"
std::string DescribeSection(const SectionMeasures& measures);
