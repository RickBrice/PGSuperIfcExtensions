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
#include "IfcProjectMappingTable.h"
#include "MenuLogo.h"

#include <IFace/Tools.h>
#include <IFace/Project.h>
#include <EAF/EAFDocument.h>
#include <EAF/EAFUIIntegration.h>
#include <EAF/AutoProgress.h>

#include <fstream>
#include <map>
#include <string>

CPGSuperDataExporter::CPGSuperDataExporter()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   LoadMenuLogo(m_Bitmap);
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

   // The Mapping Table tab defaults to the project's table (its last export table, else its import table), else the
   // configured table
   bool bFromImport = false;
   std::filesystem::path project_table = GetProjectExportTable(pBroker, &bFromImport);
   if (!project_table.empty())
   {
      CString label(bFromImport ? _T("This project's IFC import table: ") : _T("This project's last IFC export table: "));
      options_sheet.m_TablePage.SetDefault(label + GetMappingTableDisplayName(project_table), project_table, MappingTableSource::Project);
   }

   // Round trips: the IFC tab offers to update a copy of the model the project was imported from, if it's unchanged
   GET_IFACE2(pBroker, IIfcRoundTrip, pRoundTrip);
   const auto& registry = pRoundTrip->GetRoundTripRegistry();
   if (!registry.IsEmpty())
   {
      auto& page = options_sheet.m_IfcPage;
      page.m_strSourceModel = registry.source.path.c_str();
      switch (CheckSourceModel(registry.source, registry.source.path))
      {
      case SourceModelState::Unchanged: page.m_bUpdateSourceModel = true; break;
      case SourceModelState::Missing: page.m_strSourceModelProblem = _T("The model isn't there."); break;
      case SourceModelState::Changed: page.m_strSourceModelProblem = _T("The model has changed since the import."); break;
      }
   }

   if (options_sheet.DoModal() == IDCANCEL)
      return S_OK;

   GET_IFACE2(pBroker, IEAFDocument, pDoc);
   CString file_title = pDoc->GetFileTitle();
   CString file_root = pDoc->GetFileRoot();

   CString source_model = options_sheet.m_IfcPage.m_bUpdateSourceModel ? options_sheet.m_IfcPage.m_strSourceModel : CString();

   // 2. IFC file. A copy of the source model is named after it, so it isn't confused with the source
   CString ifc_default_name;
   if (source_model.IsEmpty())
   {
      ifc_default_name.Format(_T("%s%s.ifc"), file_root, file_title);
   }
   else
   {
      std::filesystem::path source(source_model.GetString());
      ifc_default_name.Format(_T("%s%s (%s).ifc"), file_root, source.stem().c_str(), file_title);
   }
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
   // a round trip uses the default options: the IFC tab's others don't apply to it, and are disabled
   CIfcExportOptions ifc_options = source_model.IsEmpty() ? options_sheet.m_IfcPage.options : CIfcExportOptions();
   ifc_options.mapping_file = options_sheet.m_TablePage.m_strTable;
   ifc_options.mapping_source = options_sheet.m_TablePage.m_Source;
   ifc_options.source_model = source_model;
   if (bExportIds)
      ifc_options.segment_global_ids = &segment_global_ids;

   CIfcExporter ifc_builder;
   bool bIfcResult = false;
   CString strError;
   try
   {
      bIfcResult = ifc_builder.BuildModel(pBroker, ifc_options, ifc_path);
   }
   catch (const std::exception& e)
   {
      // e.g. the mapping table can't be used (CIfcMappingTableException). The message says what to do
      strError = e.what();
   }

   // the project remembers the table, so its next export defaults to it
   if (bIfcResult)
   {
      try
      {
         RememberExportTable(pBroker, *CIfcMappingTable::LoadActive(std::filesystem::path(ifc_options.mapping_file.GetString()), ifc_options.mapping_source));
      }
      catch (const std::exception&)
      {
         // the export used the table, so it can be loaded. If it changed since, the project's table is unchanged
      }
   }

   CString strMsg;
   strMsg.Format(_T("Model export %s for %s"), (bIfcResult ? _T("successful") : _T("failed")), ifc_path);
   if (!strError.IsEmpty())
      strMsg += _T("\n\n") + strError;

   // a round trip says what it changed, and writes the details beside the IFC file
   if (bIfcResult && !source_model.IsEmpty())
   {
      const auto& merge = ifc_builder.GetMergeResult();
      std::filesystem::path log_path(ifc_path.GetString());
      log_path.replace_extension(".export.log");
      {
         std::ofstream log(log_path);
         log << "Updated a copy of " << PathToString(std::filesystem::path(source_model.GetString())) << std::endl;
         for (const auto& line : merge.log)
            log << line << std::endl;
      }

      CString strMerge;
      strMerge.Format(_T("\n\nThis is a copy of %s with PGSuper's design in it (the model itself isn't changed): %d elements updated (%d values changed, %d properties and classifications added) and %d elements added. The details are in %s."),
         source_model, (int)merge.nUpdated, (int)merge.nPropertiesChanged, (int)merge.nPropertiesAdded, (int)merge.nAdded, log_path.c_str());
      strMsg += strMerge;
   }

   // 5. Build the IDS from the same model, pinned to the just-written beam GlobalIds.
   if (bExportIds && bIfcResult)
   {
      pProgress->UpdateMessage(_T("Exporting IDS specification"));

      GET_IFACE2_NOCHECK(pBroker, IProjectProperties, pProjectProperties);
      CString bridge_name = pProjectProperties->GetBridgeName();
      if (bridge_name.IsEmpty()) bridge_name = _T("PGSuper Bridge");
      ids_options.title.Format(_T("%s - girder concrete strength requirements"), bridge_name);
      ids_options.global_id_by_segment = segment_global_ids;
      if (!source_model.IsEmpty())
      {
         // the girders keep the source's names, so they're pinned by their GlobalIds. What the copy has of the model's
         // (materials, types, strands, bars) isn't PGSuper's
         const auto& merge = ifc_builder.GetMergeResult();
         ids_options.beam_id = CIdsExportOptions::BeamId::GlobalId;
         ids_options.round_trip = true;
         ids_options.include_concrete_material = false;
         ids_options.include_strands = ids_options.include_strands && !merge.bModelTendons;
         ids_options.include_rebar = ids_options.include_rebar && !merge.bModelBars;
      }

      // the IDS checks the values where the IFC export put them
      ids_options.mapping_file = ifc_options.mapping_file;
      ids_options.mapping_source = ifc_options.mapping_source;

      CIdsExporter ids_builder;
      bool bIdsResult = false;
      CString strIdsError;
      try
      {
         bIdsResult = ids_builder.BuildSpecification(pBroker, ids_options, ids_path);
      }
      catch (const std::exception& e)
      {
         strIdsError = e.what(); // e.g. the mapping table can't be used
      }

      CString idsMsg;
      idsMsg.Format(_T("\nIDS export %s for %s"), (bIdsResult ? _T("successful") : _T("failed")), ids_path);
      strMsg += idsMsg;
      if (!strIdsError.IsEmpty())
         strMsg += _T("\n\n") + strIdsError;
      if (!bIdsResult) bIfcResult = false; // reflect the overall outcome in the icon
   }

   AfxMessageBox(strMsg, MB_OK | (bIfcResult ? MB_ICONEXCLAMATION : MB_ICONSTOP));
   return S_OK;
}
