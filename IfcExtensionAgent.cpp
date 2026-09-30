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

// IfcExtensionAgent.cpp : Implementation of CIfcExtensionAgent

#include "stdafx.h"
#include "IfcExtensions.h"
#include "IfcExtensionAgent.h"

#include <IFace\Tools.h>
#include <IFace\EditByUI.h>
#include <EAF\Transaction.h>
#include "EditGeoreferencing.h"
#include "IfcCommandLineInfo.h"
#include "IfcImporter.h"
#include "IfcExporter.h"
#include "IdsExporter.h"
#include "IfcTableIds.h"
#include "MappingTableDlg.h"
#include "IfcTableFormat.h"

#include <EAF\EAFApp.h>
#include <EAF\EAFDocument.h>
#include <EAF\EAFUtilities.h>

BEGIN_MESSAGE_MAP(CIfcExtensionAgent,CCmdTarget)
   ON_COMMAND(ID_EDIT_GEOREFERENCING,&CIfcExtensionAgent::OnEditGeoreferencing)
   ON_COMMAND(ID_OPTIONS_IFC_MAPPING_TABLE,&CIfcExtensionAgent::OnIfcMappingTable)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////
// IAgentEx

bool CIfcExtensionAgent::RegisterInterfaces()
{
   EAF_AGENT_REGISTER_INTERFACES;
   REGISTER_INTERFACE(IGeoreferencing);

   return true;
}

bool CIfcExtensionAgent::Init()
{
   EAF_AGENT_INIT;
   CREATE_LOGFILE(_T("IfcExtensionAgent"));

   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   VERIFY(m_bmpMenu.LoadBitmap(IDB_BSI));

   return true;
}

bool CIfcExtensionAgent::Reset()
{
   EAF_AGENT_RESET;
   m_bmpMenu.DeleteObject();
   return true;
}

bool CIfcExtensionAgent::ShutDown()
{
   EAF_AGENT_SHUTDOWN;
   return true;
}

CLSID CIfcExtensionAgent::GetCLSID() const
{
   return CLSID_PGSuperIfcExtensionAgent;
}

////////////////////////////////////////////////////////////////////
// IAgentPersist

WBFL::EAF::Broker::LoadResult CIfcExtensionAgent::Load(WBFL::System::IStructuredLoad* pStrLoad)
{
   if ( !pStrLoad->BeginUnit(_T("IfcExtensionAgent")) )
      return WBFL::EAF::Broker::LoadResult::Error;

   m_GeoreferencingData.Load(pStrLoad);

   if ( !pStrLoad->EndUnit() )
      return WBFL::EAF::Broker::LoadResult::Error;

   return WBFL::EAF::Broker::LoadResult::Success;
}

bool CIfcExtensionAgent::Save(WBFL::System::IStructuredSave* pStrSave)
{
   pStrSave->BeginUnit(_T("IfcExtensionAgent"),1.0);
   m_GeoreferencingData.Save(pStrSave);
   pStrSave->EndUnit();
   return true;
}

////////////////////////////////////////////////////////////////////
// IAgentUIIntegration

bool CIfcExtensionAgent::IntegrateWithUI(bool bIntegrate)
{
   if ( bIntegrate )
   {
      RegisterUIExtensions();
      CreateMenus();
   }
   else
   {
      RemoveMenus();
      UnregisterUIExtensions();
   }

   return true;
}

void CIfcExtensionAgent::CreateMenus()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());

   GET_IFACE(IEAFMainMenu,pMainMenu);
   auto pMenu = pMainMenu->GetMainMenu();

   UINT editPos = pMenu->FindMenuItem(_T("Edit"));
   m_pEditMenu = pMenu->GetSubMenu(editPos);

   UINT alignmentPos = m_pEditMenu->FindMenuItem(_T("Alignment..."));

   auto callback = std::dynamic_pointer_cast<WBFL::EAF::ICommandCallback>(shared_from_this());
   m_pEditMenu->InsertMenu(alignmentPos, ID_EDIT_GEOREFERENCING, _T("&Georeferencing..."), callback);
   m_pEditMenu->SetMenuItemBitmaps(ID_EDIT_GEOREFERENCING, MF_BYCOMMAND, &m_bmpMenu, nullptr, callback);

   UINT optionsPos = pMenu->FindMenuItem(_T("Options"));
   if (optionsPos != (UINT)-1)
   {
      m_pOptionsMenu = pMenu->GetSubMenu(optionsPos);
      m_pOptionsMenu->AppendMenu(ID_OPTIONS_IFC_MAPPING_TABLE, _T("&IFC Mapping Table..."), callback);
      m_pOptionsMenu->SetMenuItemBitmaps(ID_OPTIONS_IFC_MAPPING_TABLE, MF_BYCOMMAND, &m_bmpMenu, nullptr, callback);
   }
}

void CIfcExtensionAgent::RemoveMenus()
{
   if ( m_pEditMenu )
   {
      auto callback = std::dynamic_pointer_cast<WBFL::EAF::ICommandCallback>(shared_from_this());
      m_pEditMenu->RemoveMenu(ID_EDIT_GEOREFERENCING, MF_BYCOMMAND, callback);
      m_pEditMenu.reset();
   }

   if ( m_pOptionsMenu )
   {
      auto callback = std::dynamic_pointer_cast<WBFL::EAF::ICommandCallback>(shared_from_this());
      m_pOptionsMenu->RemoveMenu(ID_OPTIONS_IFC_MAPPING_TABLE, MF_BYCOMMAND, callback);
      m_pOptionsMenu.reset();
   }
}

////////////////////////////////////////////////////////////////////
// ICommandCallback

BOOL CIfcExtensionAgent::OnCommandMessage(UINT nID, int nCode, void* pExtra, AFX_CMDHANDLERINFO* pHandlerInfo)
{
   return OnCmdMsg(nID, nCode, pExtra, pHandlerInfo);
}

BOOL CIfcExtensionAgent::GetStatusBarMessageString(UINT nID, CString& rMessage) const
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());

   if (rMessage.LoadString(nID))
   {
      // first newline terminates actual string
      rMessage.Replace('\n', '\0');
   }
   else
   {
      TRACE1("Warning (CIfcExtensionAgent): no message line prompt for ID 0x%04X.\n", nID);
   }

   return TRUE;
}

BOOL CIfcExtensionAgent::GetToolTipMessageString(UINT nID, CString& rMessage) const
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   CString string;
   if (string.LoadString(nID))
   {
      // tip is after first newline
      int pos = string.Find('\n');
      if (0 < pos)
         rMessage = string.Mid(pos + 1);
   }
   else
   {
      TRACE1("Warning (CIfcExtensionAgent): no tool tip for ID 0x%04X.\n", nID);
   }

   return TRUE;
}

void CIfcExtensionAgent::OnEditGeoreferencing()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   GET_IFACE(IEditByUI,pEditByUI);
   pEditByUI->EditAlignmentDescription(_T("Georeferencing"));
}

void CIfcExtensionAgent::OnIfcMappingTable()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   CMappingTableDlg dlg(EAFGetMainFrame());
   dlg.m_strTable = CIfcMappingTable::GetTableSetting().c_str();
   if (dlg.DoModal() == IDOK)
      CIfcMappingTable::SetTableSetting(std::filesystem::path(dlg.m_strTable.GetString()));
}

void CIfcExtensionAgent::RegisterUIExtensions()
{
   // RegisterEditAlignmentCallback is declared on the shared IExtendUI base, but only
   // IExtendPGSuperUI/IExtendPGSpliceUI are ever registered on the broker (see
   // CPGSuperDocProxyAgent::RegisterInterfaces) - IID_IExtendUI itself is never
   // queryable directly. This plugin is PGSuper-specific, so IExtendPGSuperUI is the
   // right one to go through (same pattern ExampleExtensionAgent uses).
   GET_IFACE(IExtendPGSuperUI,pExtendPGSuperUI);
   m_EditAlignmentCallbackID = pExtendPGSuperUI->RegisterEditAlignmentCallback(this);
}

void CIfcExtensionAgent::UnregisterUIExtensions()
{
   GET_IFACE(IExtendPGSuperUI,pExtendPGSuperUI);
   pExtendPGSuperUI->UnregisterEditAlignmentCallback(m_EditAlignmentCallbackID);
}

////////////////////////////////////////////////////////////////////
// IGeoreferencing
void CIfcExtensionAgent::SetGeoreferencingData(const GeoreferencingData& data)
{
   m_GeoreferencingData = data;
}

const GeoreferencingData& CIfcExtensionAgent::GetGeoreferencingData() const
{
   return m_GeoreferencingData;
}


////////////////////////////////////////////////////////////////////
// IEditAlignmentCallback

CPropertyPage* CIfcExtensionAgent::CreatePropertyPage(IEditAlignmentData* pAlignmentData)
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   CGeoReferencingPage* pPage = new CGeoReferencingPage();
   pPage->m_GeoRefData = m_GeoreferencingData;
   return pPage;
}

std::unique_ptr<WBFL::EAF::Transaction> CIfcExtensionAgent::OnOK(CPropertyPage* pPage,IEditAlignmentData* pAlignmentData)
{
   CGeoReferencingPage* pMyPage = (CGeoReferencingPage*)pPage;
   auto pTxn = std::make_unique<txnEditGeoreferencing>(m_GeoreferencingData, pMyPage->m_GeoRefData);
   m_GeoreferencingData = pMyPage->m_GeoRefData;
   return pTxn;
}


////////////////////////////////////////////////////////////////////
// IEAFProcessCommandLine
BOOL CIfcExtensionAgent::ProcessCommandLineOptions(CEAFCommandLineInfo& cmdInfo)
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());

   // cmdInfo is the command line information from the application. The application
   // doesn't know about this extension at the time the command line parameters are parsed
   //
   // Re-parse the parameters with our own command line information object
   CIfcCommandLineInfo ifcCmdInfo;
   EAFGetApp()->ParseCommandLine(ifcCmdInfo);
   if (!ifcCmdInfo.m_bIfcImport && !ifcCmdInfo.m_bIfcExport && !ifcCmdInfo.m_bTableToIds && !ifcCmdInfo.m_bIdsToTable && !ifcCmdInfo.m_bFormatTable)
      return FALSE; // not our command line

   if (ifcCmdInfo.m_bError)
   {
      ifcCmdInfo.SetErrorInfo(ifcCmdInfo.GetUsageMessage());
      cmdInfo = ifcCmdInfo;
      return TRUE; // our command line, but it isn't correct. The application reports the error
   }

   cmdInfo = ifcCmdInfo; // copies m_bCommandLineMode so the application shuts down when we are done

   if (ifcCmdInfo.m_bIfcImport)
      ImportFromCommandLine(ifcCmdInfo);
   else if (ifcCmdInfo.m_bIfcExport)
      ExportFromCommandLine(ifcCmdInfo);
   else if (ifcCmdInfo.m_bTableToIds)
      TableToIdsFromCommandLine(ifcCmdInfo);
   else if (ifcCmdInfo.m_bIdsToTable)
      IdsToTableFromCommandLine(ifcCmdInfo);
   else
      FormatTableFromCommandLine(ifcCmdInfo);

   return TRUE;
}

void CIfcExtensionAgent::ImportFromCommandLine(const CIfcCommandLineInfo& ifcCmdInfo)
{
   // The template given on the command line has been opened as a new project.
   // Import the IFC model into it
   CIfcImportOptions options;
   options.interactive = false;
   options.log_file = ifcCmdInfo.m_strLogFile;
   options.mapping_file = ifcCmdInfo.m_strMappingFile;
   CString strIfcFile(ifcCmdInfo.m_strIfcFile);
   HRESULT hr = CIfcImporter(EAFGetBroker()).ImportFromIFC(strIfcFile, options);

   // Save the project even if the import failed so the partial result can be inspected
   CEAFDocument* pDoc = EAFGetDocument();
   BOOL bSaved = pDoc->DoSave(ifcCmdInfo.m_strOutFile, TRUE);

   std::wofstream log(ifcCmdInfo.m_strLogFile.GetString(), std::ios::app);
   log << (SUCCEEDED(hr) ? _T("IFC import succeeded") : _T("IFC import failed")) << std::endl;
   log << (bSaved ? _T("Project saved to ") : _T("Unable to save project to ")) << ifcCmdInfo.m_strOutFile.GetString() << std::endl;
}

void CIfcExtensionAgent::ExportFromCommandLine(const CIfcCommandLineInfo& ifcCmdInfo)
{
   // The project given on the command line has been opened.
   // Export it with the default export options
   CIfcExportOptions options;
   options.display_units_for_properties = ifcCmdInfo.m_bDisplayUnitsForProperties;
   options.mapping_file = ifcCmdInfo.m_strMappingFile;

   bool bResult = false;
   CString strError;
   CIfcExporter exporter;
   try
   {
      bResult = exporter.BuildModel(EAFGetBroker(), options, ifcCmdInfo.m_strIfcFile);
   }
   catch (const std::exception& e)
   {
      strError = e.what();
   }

   // the design-value IDS for the exported model, with the default options
   bool bIdsResult = false;
   CString strIdsError;
   if (bResult && !ifcCmdInfo.m_strIdsFile.IsEmpty())
   {
      try
      {
         CIdsExportOptions ids_options;
         ids_options.enabled = true;
         ids_options.mapping_file = ifcCmdInfo.m_strMappingFile;
         bIdsResult = CIdsExporter().BuildSpecification(EAFGetBroker(), ids_options, ifcCmdInfo.m_strIdsFile);
      }
      catch (const std::exception& e)
      {
         strIdsError = e.what();
      }
   }

   std::wofstream log(ifcCmdInfo.m_strLogFile.GetString());
   log << _T("Exported ") << ifcCmdInfo.m_strFileName.GetString() << _T(" with property values in ") << (options.display_units_for_properties ? _T("display units") : _T("system units")) << std::endl;
   for (const auto& table_file : exporter.GetMappingTableFiles())
      log << _T("IFC mapping table ") << CString(table_file.c_str()).GetString() << std::endl;
   if (!strError.IsEmpty())
      log << _T("IFC export failed: ") << strError.GetString() << std::endl;
   log << (bResult ? _T("IFC export succeeded: ") : _T("IFC export failed: ")) << ifcCmdInfo.m_strIfcFile.GetString() << std::endl;
   if (!ifcCmdInfo.m_strIdsFile.IsEmpty())
   {
      if (!strIdsError.IsEmpty())
         log << _T("IDS export failed: ") << strIdsError.GetString() << std::endl;
      log << (bIdsResult ? _T("IDS export succeeded: ") : _T("IDS export failed: ")) << ifcCmdInfo.m_strIdsFile.GetString() << std::endl;
   }
}

void CIfcExtensionAgent::TableToIdsFromCommandLine(const CIfcCommandLineInfo& ifcCmdInfo)
{
   std::ofstream log(ifcCmdInfo.m_strLogFile.GetString());
   try
   {
      std::filesystem::path table_path(ifcCmdInfo.m_strMappingFile.GetString());
      auto table = CIfcMappingTable::LoadActive(table_path);
      for (const auto& file : table->GetFiles())
         log << "IFC mapping table \"" << file.name << "\" (version " << file.version << "): " << PathToString(file.path) << ", chosen by " << MappingTableSourceDescription(file.source) << std::endl;

      std::ofstream ids(ifcCmdInfo.m_strToolIdsFile.GetString(), std::ios::binary);
      if (!ids)
         throw std::runtime_error("The IDS file can't be written: " + PathToString(std::filesystem::path(ifcCmdInfo.m_strToolIdsFile.GetString())));

      std::vector<std::string> notes;
      WriteTableAsIds(*table, ids, notes);
      for (const auto& note : notes)
         log << "Note: " << note << std::endl;
      log << "IDS written: " << PathToString(std::filesystem::path(ifcCmdInfo.m_strToolIdsFile.GetString())) << std::endl;
   }
   catch (const std::exception& e)
   {
      log << "IDS not written:" << std::endl << e.what() << std::endl;
   }
}

void CIfcExtensionAgent::FormatTableFromCommandLine(const CIfcCommandLineInfo& ifcCmdInfo)
{
   std::ofstream log(ifcCmdInfo.m_strLogFile.GetString());
   try
   {
      std::filesystem::path in_path(ifcCmdInfo.m_strFormatTableFile.GetString());
      std::filesystem::path out_path(ifcCmdInfo.m_strTableFile.GetString());

      std::ifstream in(in_path, std::ios::binary);
      if (!in)
         throw std::runtime_error("The mapping table can't be read: " + PathToString(in_path));
      std::stringstream content;
      content << in.rdbuf();
      auto table = nlohmann::ordered_json::parse(content.str()); // throws with the line and column

      std::ofstream out(out_path, std::ios::binary);
      if (!out)
         throw std::runtime_error("The mapping table can't be written: " + PathToString(out_path));
      out << FormatMappingTable(table);

      log << "Mapping table " << PathToString(in_path) << " written as " << PathToString(out_path) << std::endl;
   }
   catch (const std::exception& e)
   {
      log << "Mapping table not written:" << std::endl << e.what() << std::endl;
   }
}

void CIfcExtensionAgent::IdsToTableFromCommandLine(const CIfcCommandLineInfo& ifcCmdInfo)
{
   std::ofstream log(ifcCmdInfo.m_strLogFile.GetString());
   try
   {
      // targets come from the standard table where the IDS and the binding file don't say
      auto standard = CIfcMappingTable::Load(std::filesystem::path(), MappingTableSource::InstalledStandard);

      std::filesystem::path ids_path(ifcCmdInfo.m_strToolIdsFile.GetString());
      std::filesystem::path table_path(ifcCmdInfo.m_strTableFile.GetString());
      std::filesystem::path binding_path(ifcCmdInfo.m_strBindingFile.GetString());
      auto result = GenerateTableFromIds(ids_path, binding_path, *standard, "", ifcCmdInfo.m_bTableExtendsStandard);

      {
         std::ofstream table(table_path, std::ios::binary);
         if (!table)
            throw std::runtime_error("The mapping table can't be written: " + PathToString(table_path));
         table << result.table_json << std::endl;
      }

      log << "Mapping table generated from " << PathToString(ids_path) << (binding_path.empty() ? std::string() : " with the binding file " + PathToString(binding_path)) << ": " << PathToString(table_path) << std::endl;
      log << result.bound << " property facets with a target, " << result.unbound << " without" << std::endl;
      for (const auto& line : result.report)
         log << line << std::endl;

      // the generated table must be one the importer and exporter accept
      try
      {
         CIfcMappingTable::Load(table_path, MappingTableSource::CommandLine);
         log << "The generated table is valid." << std::endl;
      }
      catch (const CIfcMappingTableException& e)
      {
         log << "The generated table isn't valid:" << std::endl << e.what() << std::endl;
      }
   }
   catch (const std::exception& e)
   {
      log << "Mapping table not generated:" << std::endl << e.what() << std::endl;
   }
}
