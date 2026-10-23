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

#include "resource.h"
#include <CommCtrl.h>
#pragma comment(lib, "comctl32.lib") // LoadIconWithScaleDown

// The bSI logo for a menu item, at the menu check-mark size (SM_CXMENUCHECK: 13 px on Windows 10 and 15 px on
// Windows 11 at 96 DPI, larger with display scaling), so Windows doesn't stretch it. At 13 and 15 px it's pixel art
// (IDB_BSI, IDB_BSI_15: one link color per pixel; the anti-aliased icon blends the links into mud at these sizes);
// at other sizes it's the multi-size icon (IDI_BSI) drawn on white. Call in the extension's module state
inline void LoadMenuLogo(CBitmap& bitmap)
{
   int size = ::GetSystemMetrics(SM_CXMENUCHECK);
   if (size == 13 || size == 15)
   {
      VERIFY(bitmap.LoadBitmap(size == 13 ? IDB_BSI : IDB_BSI_15));
      return;
   }

   HICON hIcon = nullptr;
   if (FAILED(::LoadIconWithScaleDown(AfxFindResourceHandle(MAKEINTRESOURCE(IDI_BSI), RT_GROUP_ICON), MAKEINTRESOURCE(IDI_BSI), size, size, &hIcon)))
   {
      VERIFY(bitmap.LoadBitmap(IDB_BSI));
      return;
   }

   CWindowDC screen(nullptr);
   CDC dc;
   dc.CreateCompatibleDC(&screen);
   bitmap.CreateCompatibleBitmap(&screen, size, size);
   CBitmap* pOld = dc.SelectObject(&bitmap);
   dc.FillSolidRect(0, 0, size, size, RGB(255, 255, 255));
   ::DrawIconEx(dc, 0, 0, hIcon, size, size, 0, nullptr, DI_NORMAL);
   dc.SelectObject(pOld);
   ::DestroyIcon(hIcon);
}
