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

#include "Geometry.h"
#include <map>
#include <set>
#include <memory>

namespace WBFL { namespace EAF { class Broker; }; };

// The plan center (of the bounding box) of each element whose geometry can be processed, by id
std::map<int, Eigen::Vector2d> get_plan_centers(ifcopenshell::file& file, const std::set<int>& ids);

// Picks the structure to import when the bridge is several structures side by side: deck slabs whose offset ranges
// don't overlap, with no slab across the joint between them (e.g. parallel bridges for the two directions of a highway
// in one IfcBridge, each with its own deck, girders, bearings, and barriers). PGSuper models one structure, so the one
// on the alignment (its deck spans offset 0) is imported. Each element of an element role goes with the structure whose
// deck is over its plan center (else the structure nearest in offset), and the elements of the other structures are
// left out of the roles for the rest of the import (CIfcImporter::ExcludeElements). The decision is logged.
// Nothing is left out when the deck is one structure, or when no structure is on the alignment.
// Call it after the alignment is known and before the elements are used
void select_structure(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file);
