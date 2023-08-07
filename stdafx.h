///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright © 1999-2023  Washington State Department of Transportation
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

// stdafx.h : include file for standard system include files,
//      or project specific include files that are used frequently,
//      but are changed infrequently

#pragma once

#include <PGSuperAll.h>
#include "resource.h"

#include <atlbase.h>
//You may derive a class from CComModule and use it if you want to override
//something, but do not change the name of _Module
extern CComModule _Module;
#include <atlcom.h>

#include <AgentTools.h>

#pragma warning(disable:4250) // 
#include <ifcparse/IfcHierarchyHelper.h>
//#include <ifcparse/Ifc4x1.h>
//#include <ifcparse/Ifc4x2.h>
//#include <ifcparse/Ifc4x3_rc1.h>
//#include <ifcparse/Ifc4x3_rc2.h>
//#include <ifcparse/Ifc4x3_rc3.h>
//#include <ifcparse/Ifc4x3_rc4.h>
#include <ifcparse/Ifc4x3_tc1.h>
#include <ifcparse/Ifc4x3_add1.h>
