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

#include "stdafx.h"
#include "EditGeoreferencing.h"

#include <IFace\Tools.h>
#include <EAF\EAFUtilities.h>

txnEditGeoreferencing::txnEditGeoreferencing(const GeoreferencingData& oldData, const GeoreferencingData& newData)
{
   m_Data[0] = oldData;
   m_Data[1] = newData;
}

bool txnEditGeoreferencing::Execute()
{
   Execute(1);
   return true;
}

void txnEditGeoreferencing::Undo()
{
   Execute(0);
}

void txnEditGeoreferencing::Execute(int i)
{
   auto pBroker = EAFGetBroker();
   GET_IFACE2(pBroker,IGeoreferencing,pGeoref);
   pGeoref->SetGeoreferencingData(m_Data[i]);
}

std::unique_ptr<WBFL::EAF::Transaction> txnEditGeoreferencing::CreateClone() const
{
   return std::make_unique<txnEditGeoreferencing>(m_Data[0], m_Data[1]);
}

std::_tstring txnEditGeoreferencing::Name() const
{
   return _T("Edit Georeferencing");
}

bool txnEditGeoreferencing::IsUndoable() const
{
   return true;
}

bool txnEditGeoreferencing::IsRepeatable() const
{
   return false;
}
