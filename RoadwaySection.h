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

namespace WBFL { namespace EAF { class Broker; }; };

#include "Geometry.h"
#include "TrafficBarriers.h"

// Imports the roadway section (cross slopes) from the top of the deck (devdocs/IfcImportPlan.md, Phase 4). The deck is
// cut normal to the alignment at the quarter points of each span and near the ends of the bridge. The top of the deck
// between the barriers (the roadway) is reduced to its corners, which are the ridge points, with a ridge point at the
// alignment, where the profile grade is. The sections become the roadway section templates: one if they're all the
// same, otherwise where the section changes, extended linearly to the end piers in a transition. A section missing a
// ridge point the others have gets it; stations where the section can't be measured, or still has a different number
// of ridge points than most, are left out
void import_roadway_section(std::shared_ptr<WBFL::EAF::Broker> pBroker, const Mesh& deck, const std::vector<ElementMesh>& barriers, const std::vector<Float64>& pier_stations);

// Measures the profile from the top of the deck at the alignment (devdocs/IfcImportPlan.md, Phase 4), about every 10 ft
// between the end piers, and fits tangents and symmetric parabolic vertical curves to it within 1/8 in. With bReplace
// (the alignment has no vertical alignment) the fit becomes the profile; otherwise the deck is compared with the
// alignment's profile and a difference is logged
void import_profile(std::shared_ptr<WBFL::EAF::Broker> pBroker, const Mesh& deck, const std::vector<ElementMesh>& barriers, const std::vector<Float64>& pier_stations, bool bReplace);
