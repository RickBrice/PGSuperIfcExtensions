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
#include <array>
#include <optional>
#include <string>
#include <vector>

struct Mesh;
class GirderLibraryEntry;
interface IShape;
namespace PGS { namespace Beams { class BeamFactory; } }

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
   std::vector<double> area_profile;  // area of the section above each AREA_PROFILE_STEP below the top (exact, so it
                                      // changes smoothly with the dimensions, for fitting)
   SectionOutline outline;            // the outline measured (section coordinates)
};

constexpr double WIDTH_PROFILE_STEP = 0.003175; // 1/8 in
constexpr double AREA_PROFILE_STEP = 0.0127;    // 1/2 in

// The tolerances for a section to match another. The same for every girder family for now
// (devdocs/GirderImportPlan.md, decision 2), but kept as data so they can differ by family later
struct SectionTolerances
{
   double length = 0.003175;      // 1/8 in: height, centroid, and widths at each depth (and the depth where a width changes)
   double ratio = 0.01;           // area and moment of inertia
   double local_depth = 0.0127;   // 1/2 in: the total depth over which a fitted section's widths may differ (FittedSectionComparison)
};

// The closed loops of a mesh cut by the vertical plane through a point, normal to a horizontal axis, largest (by area)
// first. A loop's points are (distance to the left of the point, looking along the axis; elevation)
std::vector<SectionOutline::Loop> CutMesh(const Mesh& mesh, const Eigen::Vector2d& point, Eigen::Vector2d axis);

// The segments where the faces of a mesh are cut by the same plane, in the same coordinates, without joining them into
// loops: a mesh of several solids that touch (e.g. a deck slab and its haunches) can't always be joined correctly
std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> CutMeshSegments(const Mesh& mesh, const Eigen::Vector2d& point, Eigen::Vector2d axis);

// A girder section in section coordinates, but not yet upright or at the origin (e.g. as cut or projected from the
// model), turned upright and moved so the origin is at its top center. A girder can be rotated about its axis (e.g. to
// the cross slope of the deck): the section is turned, within 15 degrees, to where it's the least high
SectionOutline UprightSection(SectionOutline outline);

// The section of a girder mesh cut by a vertical plane normal to the girder axis (start to end, in plan) at a
// fraction of its length. Returns nullopt if the plane doesn't cut a closed outline
std::optional<SectionOutline> CutGirderSection(const Mesh& mesh, const Eigen::Vector2d& start, const Eigen::Vector2d& end, double fraction);

// The outline of a WBFL shape (e.g. a girder library entry's section), moved to section coordinates
SectionOutline GetShapeOutline(IShape* pShape);

// The outline of a girder library entry's section, or nullopt if its beam factory can't make one
std::optional<SectionOutline> GetLibraryEntryOutline(const GirderLibraryEntry* pEntry);

SectionMeasures MeasureSection(const SectionOutline& outline);

// e.g. "height 72.000 in, area 913.4 in^2, I 545,894 in^4, top width 49.000 in"
std::string DescribeSection(const SectionMeasures& measures);

// How a section is compared with another (devdocs/GirderImportPlan.md, G1). Comparisons are strategies, so other
// approaches can be tried without changing the import
class ISectionComparison
{
public:
   virtual ~ISectionComparison() = default;

   // The differences the comparison doesn't allow, as text (e.g. "height 36.001 in vs 35.000 in"). Empty if the sections agree
   virtual std::vector<std::string> Compare(const SectionMeasures& model, const SectionMeasures& other) const = 0;

   // The differences the comparison allows but that are worth reporting (e.g. the corners of a flange), as text
   virtual std::vector<std::string> Allowed(const SectionMeasures& model, const SectionMeasures& other) const { return {}; }

   // A number that grows smoothly with the difference, about 1 at the tolerances: it ranks sections and is the
   // objective of fitting
   virtual double Difference(const SectionMeasures& model, const SectionMeasures& other) const = 0;
};

// The sections are the same: height, centroid, and the width at every depth within the length tolerance (a width may
// come from a depth within the tolerance, so a flange thickness that differs by less isn't a difference), area and I
// within the ratio. For recognizing an existing library entry
class ExactSectionComparison : public ISectionComparison
{
public:
   ExactSectionComparison(const SectionTolerances& tolerances = SectionTolerances()) : m_Tolerances(tolerances) {}
   std::vector<std::string> Compare(const SectionMeasures& model, const SectionMeasures& other) const override;
   double Difference(const SectionMeasures& model, const SectionMeasures& other) const override;

protected:
   SectionTolerances m_Tolerances;
};

// A fitted section represents the model's: as ExactSectionComparison, except that the widths may differ over a total
// depth up to SectionTolerances::local_depth, or over any depth if the outlines are within the length tolerance of each
// other (OutlineDistance). These are local details a girder family can't make (e.g. the corners of a flange tip), and
// small offsets on shallow tapers (e.g. a model's flange tips a tenth of an inch lower on one side, from a crown),
// where a small vertical offset is a large width difference; they're reported by Allowed. For accepting a library entry with fitted dimensions.
// Its Difference, the fitting objective, is the area between the two outlines (their symmetric difference), with the
// height, centroid, area, and I: it changes continuously with every dimension, including flange tapers and fillets
class FittedSectionComparison : public ExactSectionComparison
{
public:
   using ExactSectionComparison::ExactSectionComparison;
   std::vector<std::string> Compare(const SectionMeasures& model, const SectionMeasures& other) const override;
   std::vector<std::string> Allowed(const SectionMeasures& model, const SectionMeasures& other) const override;
   double Difference(const SectionMeasures& model, const SectionMeasures& other) const override;
};

// The area between two outlines (the area of their symmetric difference)
double OutlineDifferenceArea(const SectionOutline& a, const SectionOutline& b);

// True if the outline's coordinates are finite and none of its loops crosses itself
bool IsSimpleOutline(const SectionOutline& outline);

// The corners of an open polyline: points closer than the tolerance to the line between their neighboring corners are
// dropped (Douglas-Peucker). The ends are kept
std::vector<Eigen::Vector2d> SimplifyPolyline(const std::vector<Eigen::Vector2d>& points, double tolerance);

// The corners of a closed loop: points closer than the tolerance to the line between their neighboring corners are
// dropped (Douglas-Peucker). The first point is kept
SectionOutline::Loop SimplifyLoop(const SectionOutline::Loop& loop, double tolerance);

// The largest distance from a point of either outline to the other outline (the Hausdorff distance)
double OutlineDistance(const SectionOutline& a, const SectionOutline& b);

using GirderDimensions = std::vector<std::pair<std::_tstring, Float64>>;

// The dimensions of a girder family that give a section closest to the model's section (devdocs/GirderImportPlan.md,
// G1 step 3)
struct DimensionFit
{
   GirderDimensions dimensions;
   SectionMeasures section;              // the section of the fitted dimensions
   double difference = DBL_MAX;          // ISectionComparison::Difference
   std::vector<std::string> differences; // ISectionComparison::Compare: empty if the fitted section is accepted
   std::vector<std::string> allowed;     // ISectionComparison::Allowed
};

// Fits the dimensions of a library entry's beam factory to a section. Fitting is a strategy, like ISectionComparison
class IDimensionFitter
{
public:
   virtual ~IDimensionFitter() = default;

   // The fitted dimensions of a beam factory, starting from the given dimensions (e.g. a library entry's, or the
   // factory's defaults), judged by the comparison. nullopt if the factory can't make a section
   virtual std::optional<DimensionFit> Fit(std::shared_ptr<PGS::Beams::BeamFactory> factory, const GirderDimensions& start, const SectionMeasures& model, const ISectionComparison& comparison) const = 0;
};

// The fitters share these rules: the length dimensions are fitted (scalars and booleans keep their values), dimensions
// that don't change the section (e.g. end blocks) keep their values, dimensions the factory doesn't accept, with a
// negative length, or whose section isn't simple (IsSimpleOutline: a factory may accept dimensions that make an outline
// cross itself) are rejected, and a fit has a budget (evaluations and time)

// Only the height: the dimension that changes the height the most (e.g. H, or the web height) is solved for the model's
// height. For a library entry whose flanges match the model's (FlangesMatch): girders of a series (e.g. the WSDOT WF
// girders) have the same flanges and differ only in height
class HeightDimensionFitter : public IDimensionFitter
{
public:
   std::optional<DimensionFit> Fit(std::shared_ptr<PGS::Beams::BeamFactory> factory, const GirderDimensions& start, const SectionMeasures& model, const ISectionComparison& comparison) const override;
};

// All of the dimensions by least squares (Levenberg-Marquardt), after solving the height. The residuals are the
// average widths in each band of depth (from the exact area profile, so they change continuously with the
// dimensions), the height, the centroid, the area, and I, each about 1 at the tolerances
class LevenbergMarquardtDimensionFitter : public IDimensionFitter
{
public:
   LevenbergMarquardtDimensionFitter(int max_iterations = 100) : m_MaxIterations(max_iterations) {}
   std::optional<DimensionFit> Fit(std::shared_ptr<PGS::Beams::BeamFactory> factory, const GirderDimensions& start, const SectionMeasures& model, const ISectionComparison& comparison) const override;

private:
   int m_MaxIterations;
};

// All of the dimensions by least squares (Levenberg-Marquardt) of the distance between the outlines, after solving the
// height. Girder sections are symmetric about their vertical axis, so the left halves are compared. Two fits, and the
// better is kept: the corners of the outlines paired in order along them, then the outline distance (for polygonal
// sections, e.g. WF girders, whose corners fix the dimensions), and the outline distance alone (for curved sections,
// e.g. NU girders' fillets). The outline distance measures points spaced evenly along each outline to the nearest edge
// of the other, so a fit from a factory's defaults isn't held where a flange edge lines up with the wrong edge, as it
// is when widths are compared depth by depth. Least squares steps stop at zero lengths
class OutlineDimensionFitter : public IDimensionFitter
{
public:
   OutlineDimensionFitter(int max_iterations = 100) : m_MaxIterations(max_iterations) {}
   std::optional<DimensionFit> Fit(std::shared_ptr<PGS::Beams::BeamFactory> factory, const GirderDimensions& start, const SectionMeasures& model, const ISectionComparison& comparison) const override;

private:
   int m_MaxIterations;
};

// All of the dimensions by a downhill simplex search (restarted), then a pattern search one dimension at a time, after
// solving the height. Minimizes ISectionComparison::Difference. Kept to compare with the other fitters
class SimplexDimensionFitter : public IDimensionFitter
{
public:
   std::optional<DimensionFit> Fit(std::shared_ptr<PGS::Beams::BeamFactory> factory, const GirderDimensions& start, const SectionMeasures& model, const ISectionComparison& comparison) const override;
};

// True if the top and bottom flanges of two sections match: the widths in the top fraction of the shorter section
// measured down from the tops, and in the bottom fraction measured up from the bottoms (the top of one isn't compared
// with the bottom of the other). Girders of a series match, whatever their heights
bool FlangesMatch(const SectionMeasures& model, const SectionMeasures& other, const SectionTolerances& tolerances, double fraction = 0.4);

// The fit with a girder of constant depth. The section is fitted at one place, so dimensions that only shape the rest of
// the girder (e.g. a tapered I-beam's depth at its far end) keep their starting values, which can make the girder vary
// in depth. When the fitted girder is higher at one end than the other (BeamFactory::GetBeamHeight), a dimension that
// doesn't change the fitted section is solved so the heights are the same. nullopt if no dimension does that: the
// factory can't make the model's girder of constant depth with this section. (BeamFactory::IsPrismatic isn't used: it
// says whether the section properties are constant, so end blocks make a girder non-prismatic)
std::optional<DimensionFit> MakeConstantDepth(std::shared_ptr<PGS::Beams::BeamFactory> factory, const DimensionFit& fit);

// True if the factory's default girder has the same height at both ends (a factory made for girders of constant
// depth, unlike e.g. a tapered I-beam)
bool HasConstantDepthDefaults(const PGS::Beams::BeamFactory& factory);

// True if two sections measured from a model are the same girder's: outlines within the length tolerance of each other
// (OutlineDistance), and area, moment of inertia, and centroid within the tolerances. The girders of a span (and often
// of the whole bridge) are one shape; sections cut from their geometry differ only by the measurement (e.g. a girder
// sitting 0.1 in lower on one side under a crowned deck), so a girder with the section of one already matched uses its
// library entry without being fitted again
bool SameGirderSection(const SectionMeasures& a, const SectionMeasures& b, const SectionTolerances& tolerances = SectionTolerances());

// Generates the strand grid of a new library entry from its section (devdocs/GirderImportPlan.md, decision 5): positions
// on a 2 in grid (pairs at +/-1, 3, 5, ... in) with at least 2 in from the strand center to the section's edges and voids.
// Straight strands are in the lower half of the section; harped strands are in the innermost column that runs through
// mid-height (the web), at the top of it at the ends and at the bottom of it at the harp points; temporary strands are in
// the top quarter. Straight and temporary strands stay out of the harped column. None can be debonded. A section without
// such a column (e.g. a slab) has no harped strands. Returns the number of straight, harped, and temporary positions
std::array<size_t, 3> GenerateStrandGrid(GirderLibraryEntry& entry, const SectionOutline& outline);
