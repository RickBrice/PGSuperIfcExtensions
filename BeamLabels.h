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

// Shared derivation of the identifiers the IFC exporter writes onto girder IfcBeam
// occurrences. The IDS exporter reuses these so a project IDS can pin an <applicability>
// facet to exactly the beam the IFC exporter produced.
//
// GetIfcBeamName() below reproduces, byte for byte, the Name strings built inline in
// IfcExporter.cpp CreateBridge<Schema>() / CreateGirder<Schema>() (see the line
// references on the function). Those call sites should be switched to call this helper
// so the two exporters can never drift; until then, keep the two in sync by hand.

#include <string>
#include <sstream>
#include <atlconv.h>
#include <PsgLib\Keys.h>
#include <psgLib\GirderLabel.h>

// Name assigned to a girder segment's IfcBeam occurrence.
//
//   non-spliced girder (one segment)  -> pgsGirderLabel::GetGirderLabel(girderKey)
//                                        e.g. "Span 1, Girder B"  (unique per bridge)
//   spliced girder segment            -> "Segment <n>"  (n == 1-based segment index;
//                                        NOT unique across girders - use GlobalId to
//                                        identify spliced segments in an IDS)
//
// IMPORTANT: this must NOT force pgsGirderLabel::UseAlphaLabel() one way or the other.
// GIRDER_LABEL (IfcExporter.cpp CreateBridge<Schema>() / CreateGirder<Schema>()) uses
// whatever the ambient alpha/numeric setting is at export time - e.g. "Span 1, Girder B"
// when alpha labelling is on, "Span 1, Girder 2" when it is off - and that ambient
// setting is a user preference, not something derivable from the .pgs file. Forcing it
// here would silently desynchronize the Name this helper predicts from the Name the IFC
// exporter actually wrote, breaking every IDS <applicability> that pins on it.
inline std::string GetIfcBeamName(const CSegmentKey& segmentKey, bool bSpliced)
{
   USES_CONVERSION;

   if (bSpliced)
   {
      std::ostringstream os;
      os << "Segment " << (segmentKey.segmentIndex + 1);
      return os.str();
   }

   std::_tstring strLabel = pgsGirderLabel::GetGirderLabel(CGirderKey(segmentKey.groupIndex, segmentKey.girderIndex));
   return std::string(T2A(strLabel.c_str()));
}
