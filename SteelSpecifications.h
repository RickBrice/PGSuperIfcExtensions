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

// Material specifications of reinforcing bars, for the exported material properties (usBrPset_ACI_ReinforcingMaterial)

#include <string>

inline std::string GetRebarSpecification(const WBFL::Materials::Rebar* pRebar)
{
   std::string spec;
   switch (pRebar->GetType())
   {
   case WBFL::Materials::Rebar::Type::A615: spec = "ASTM A615 (AASHTO M31)"; break;
   case WBFL::Materials::Rebar::Type::A706: spec = "ASTM A706"; break;
   case WBFL::Materials::Rebar::Type::A1035: spec = "ASTM A1035"; break;
   default: spec = "Unknown"; break;
   }
   return spec;
}

inline std::string GetRebarSpecificationEdition(const WBFL::Materials::Rebar* pRebar)
{
   std::string edition;
   switch (pRebar->GetType())
   {
   case WBFL::Materials::Rebar::Type::A615: edition = "2026"; break;
   case WBFL::Materials::Rebar::Type::A706: edition = "2026"; break;
   case WBFL::Materials::Rebar::Type::A1035: edition = "2024"; break;
   default: edition = "Unknown"; break;
   }
   return edition;
}
