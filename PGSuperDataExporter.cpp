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

// PGSuperDataExporter.cpp : Implementation of CPGSuperDataExporter
#include "stdafx.h"
#include "IfcExtensions.h"
#include "PGSuperDataExporter.h"
#include "IfcExporter.h"
#include "IdsExporter.h"
#include "ExportOptionsSheet.h"

#include <IFace/Tools.h>
#include <IFace/Project.h>
#include <EAF/EAFDocument.h>
#include <EAF/EAFUIIntegration.h>
#include <EAF/AutoProgress.h>

#include <map>
#include <string>

CPGSuperDataExporter::CPGSuperDataExporter()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   VERIFY(m_Bitmap.LoadBitmap(IDB_BSI));
}

STDMETHODIMP CPGSuperDataExporter::Init(UINT nCmdID)
{
   return S_OK;
}

CString CPGSuperDataExporter::GetMenuText() const
{
   return CString("Bridge Model to IFC");
}

HBITMAP CPGSuperDataExporter::GetBitmapHandle() const
{
   return m_Bitmap;
}

CString CPGSuperDataExporter::GetCommandHintText() const
{
   return CString("Status line hint text\nTool tip text");
}

STDMETHODIMP CPGSuperDataExporter::Export(std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());

   // 1. Export options: one dialog, two tabs (IFC, IDS). The IDS tab's own "Enable IDS
   // export" checkbox (options.enabled) is the only switch for step 3/4/6 below.
   CExportOptionsSheet options_sheet(_T("Export Options"));
   if (options_sheet.DoModal() == IDCANCEL)
      return S_OK;

   GET_IFACE2(pBroker, IEAFDocument, pDoc);
   CString file_title = pDoc->GetFileTitle();
   CString file_root = pDoc->GetFileRoot();

   // 2. IFC file
   CString ifc_default_name;
   ifc_default_name.Format(_T("%s%s.ifc"), file_root, file_title);
   CFileDialog ifc_dlg(FALSE, _T("ifc"), ifc_default_name, OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, _T("IFC File (*.ifc)|*.ifc||"));
   if (ifc_dlg.DoModal() != IDOK)
      return S_OK;
   CString ifc_path = ifc_dlg.GetPathName();

   // 3. Project IDS file (only when enabled on the IDS tab). Collected before the build
   // so the whole flow is front-loaded.
   bool bExportIds = options_sheet.m_IdsPage.options.enabled;
   CIdsExportOptions ids_options = options_sheet.m_IdsPage.options;
   CString ids_path;
   if (bExportIds)
   {
      CString ids_default_name;
      ids_default_name.Format(_T("%s%s.ids"), file_root, file_title);
      CFileDialog ids_dlg_file(FALSE, _T("ids"), ids_default_name, OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, _T("IDS File (*.ids)|*.ids||"));
      if (ids_dlg_file.DoModal() == IDOK)
         ids_path = ids_dlg_file.GetPathName();
      else
         bExportIds = false; // skip the IDS, still export the IFC
   }

   // 4. Build the IFC model, capturing the beam GlobalIds when an IDS will follow.
   //
   // One progress window spans both builds. CIfcExporter::BuildModel() opens its own
   // (via WBFL::EAF::AutoProgress) for the IFC phase; IEAFProgress is ref-counted
   // (CEAFDocProxyAgent::CreateProgressWindow/DestroyProgressWindow), so nesting it
   // inside this outer one just increments/decrements that count instead of tearing
   // the window down and losing it before the IDS phase starts.
   GET_IFACE2(pBroker, IEAFProgress, pProgress);
   WBFL::EAF::AutoProgress ap(pProgress);
   pProgress->UpdateMessage(_T("Exporting bridge model"));

   std::map<CSegmentKey, std::string> segment_global_ids;
   CIfcExportOptions ifc_options = options_sheet.m_IfcPage.options;
   if (bExportIds)
      ifc_options.segment_global_ids = &segment_global_ids;

   CIfcExporter ifc_builder;
   bool bIfcResult = ifc_builder.BuildModel(pBroker, ifc_options, ifc_path);

   CString strMsg;
   strMsg.Format(_T("Model export %s for %s"), (bIfcResult ? _T("successful") : _T("failed")), ifc_path);

   // 5. Build the IDS from the same model, pinned to the just-written beam GlobalIds.
   if (bExportIds && bIfcResult)
   {
      pProgress->UpdateMessage(_T("Exporting IDS specification"));

      GET_IFACE2_NOCHECK(pBroker, IProjectProperties, pProjectProperties);
      CString bridge_name = pProjectProperties->GetBridgeName();
      if (bridge_name.IsEmpty()) bridge_name = _T("PGSuper Bridge");
      ids_options.title.Format(_T("%s - girder concrete strength requirements"), bridge_name);
      ids_options.global_id_by_segment = segment_global_ids;

      CIdsExporter ids_builder;
      bool bIdsResult = ids_builder.BuildSpecification(pBroker, ids_options, ids_path);

      CString idsMsg;
      idsMsg.Format(_T("\nIDS export %s for %s"), (bIdsResult ? _T("successful") : _T("failed")), ids_path);
      strMsg += idsMsg;
      if (!bIdsResult) bIfcResult = false; // reflect the overall outcome in the icon
   }

   AfxMessageBox(strMsg, MB_OK | (bIfcResult ? MB_ICONEXCLAMATION : MB_ICONSTOP));
   return S_OK;
}
