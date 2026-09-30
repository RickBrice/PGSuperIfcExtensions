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

// MappingEditorDoc.cpp : the mapping table editor's application plugin and document
//

#include "stdafx.h"
#include "resource.h"
#include "MappingEditorDoc.h"
#include "MappingEditorViews.h"
#include "MappingEditorPanes.h"
#include "IfcTableFormat.h"

#include <EAF\EAFApp.h>
#include <EAF\EAFMainFrame.h>
#include <EAF\EAFResources.h>
#include <EAF\EAFUtilities.h>

#include <fstream>
#include <sstream>

CString Utf8ToCString(const std::string& text)
{
   return CString(CA2W(text.c_str(), CP_UTF8));
}

std::string CStringToUtf8(const CString& text)
{
   return std::string(CW2A(text, CP_UTF8));
}

/////////////////////////////////////////////////////////////////////////////
// CMappingEditorStatusBar

void CMappingEditorStatusBar::GetStatusIndicators(const UINT** lppIDArray, int* pnIDCount)
{
   static UINT indicators[] =
   {
      ID_SEPARATOR, // status line
      EAFID_INDICATOR_AUTOSAVE_ON,
      EAFID_INDICATOR_MODIFIED,
      ID_INDICATOR_CAPS,
      ID_INDICATOR_NUM,
      ID_INDICATOR_SCRL,
   };
   *lppIDArray = indicators;
   *pnIDCount = sizeof(indicators) / sizeof(UINT);
}

/////////////////////////////////////////////////////////////////////////////
// CMappingEditorPlugin

BOOL CMappingEditorPlugin::Init(CEAFApp* pParent)
{
   return TRUE;
}

void CMappingEditorPlugin::Terminate()
{
}

void CMappingEditorPlugin::IntegrateWithUI(BOOL bIntegrate)
{
}

std::vector<CEAFDocTemplate*> CMappingEditorPlugin::CreateDocTemplates()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());

   auto pDocTemplate = new CMappingEditorDocTemplate(IDR_MAPPING_EDITOR, nullptr,
      RUNTIME_CLASS(CMappingEditorDoc), RUNTIME_CLASS(CMappingEditorFrame), RUNTIME_CLASS(CMappingTreeView), nullptr, 1);
   pDocTemplate->SetPluginApp(std::dynamic_pointer_cast<WBFL::EAF::IPluginApp>(shared_from_this()));

   return { pDocTemplate };
}

HMENU CMappingEditorPlugin::GetSharedMenuHandle()
{
   return nullptr;
}

CString CMappingEditorPlugin::GetName()
{
   return _T("IFC Mapping Table Editor");
}

CString CMappingEditorPlugin::GetDocumentationSetName()
{
   return GetName();
}

CString CMappingEditorPlugin::GetDocumentationURL()
{
   return CString(); // no online documentation yet
}

CString CMappingEditorPlugin::GetDocumentationMapFile()
{
   return CString();
}

void CMappingEditorPlugin::LoadDocumentationMap()
{
}

std::pair<WBFL::EAF::HelpResult, CString> CMappingEditorPlugin::GetDocumentLocation(LPCTSTR lpszDocSetName, UINT nID)
{
   return { WBFL::EAF::HelpResult::DocSetNotFound, CString() };
}

/////////////////////////////////////////////////////////////////////////////
// CMappingEditorDocTemplate

IMPLEMENT_DYNAMIC(CMappingEditorDocTemplate, CEAFDocTemplate)

CMappingEditorDocTemplate::CMappingEditorDocTemplate(UINT nIDResource, std::shared_ptr<WBFL::EAF::ICommandCallback> pCallback, CRuntimeClass* pDocClass,
   CRuntimeClass* pFrameClass, CRuntimeClass* pViewClass, HMENU hSharedMenu, int maxViewCount)
   : CEAFDocTemplate(nIDResource, pCallback, pDocClass, pFrameClass, pViewClass, hSharedMenu, maxViewCount)
{
   CString strDocName;
   GetDocString(strDocName, CDocTemplate::docName);

   HICON hIcon = AfxGetApp()->LoadIcon(IDR_MAPPING_EDITOR);
   m_TemplateGroup.AddItem(new CEAFTemplateItem(this, strDocName, nullptr, hIcon));
   m_TemplateGroup.SetIcon(hIcon);
}

CString CMappingEditorDocTemplate::GetTemplateGroupItemDescription(const CEAFTemplateItem* pItem) const
{
   return _T("Create an IFC mapping table: where IFC imports and exports read and write each value.");
}

/////////////////////////////////////////////////////////////////////////////
// CMappingEditorDoc

IMPLEMENT_DYNCREATE(CMappingEditorDoc, CEAFDocument)

BEGIN_MESSAGE_MAP(CMappingEditorDoc, CEAFDocument)
   ON_COMMAND(ID_MAPPING_VALIDATE, &CMappingEditorDoc::OnValidate)
END_MESSAGE_MAP()

CMappingEditorDoc::CMappingEditorDoc()
{
   EnableUIHints(FALSE);
}

void CMappingEditorDoc::TableChanged(bool bStructure)
{
   SetModifiedFlag();
   UpdateAllViews(nullptr, bStructure ? HINT_STRUCTURE : HINT_CONTENT);
}

void CMappingEditorDoc::SelectNode(const MappingNode& node)
{
   CMappingSelectHint hint;
   hint.node = node;
   UpdateAllViews(nullptr, HINT_SELECT, &hint);
}

std::filesystem::path CMappingEditorDoc::GetTablePath() const
{
   CString strPath = GetPathName();
   if (strPath.IsEmpty())
   {
      std::error_code ec;
      return std::filesystem::current_path(ec) / L"Untitled.json";
   }
   return std::filesystem::path(strPath.GetString());
}

void CMappingEditorDoc::LoadBaseTable()
{
   m_pBase.reset();
   m_BaseError.clear();

   if (!m_Table.contains("extends") || !m_Table["extends"].is_string())
      return;

   auto extends = m_Table["extends"].get<std::string>();
   if (extends.empty())
      return;

   try
   {
      if (extends == "standard")
         m_pBase = CIfcMappingTable::Load(std::filesystem::path(), MappingTableSource::InstalledStandard);
      else
         m_pBase = CIfcMappingTable::Load(GetTablePath().parent_path() / Utf8ToCString(extends).GetString(), MappingTableSource::Extends);
   }
   catch (const std::exception& e)
   {
      m_BaseError = e.what();
   }
}

bool CMappingEditorDoc::Validate(std::string& message) const
{
   std::string text = FormatMappingTable(m_Table);
   try
   {
      auto pTable = CIfcMappingTable::Load(GetTablePath(), MappingTableSource::Editor, &text);

      std::ostringstream os;
      os << "The table is valid. IFC imports and exports can use it." << std::endl << std::endl;
      bool bFirst = true;
      for (const auto& file : pTable->GetFiles())
      {
         os << (bFirst ? "" : "Extends ") << "\"" << file.name << "\" (version " << file.version << ")" << std::endl;
         os << "   " << (bFirst ? "this table" : PathToString(file.path)) << std::endl;
         bFirst = false;
      }
      message = os.str();
      return true;
   }
   catch (const std::exception& e)
   {
      message = e.what();
      return false;
   }
}

void CMappingEditorDoc::CommitPendingEdit()
{
   if (CWnd* pFrame = EAFGetMainFrame())
      pFrame->SetFocus();
}

BOOL CMappingEditorDoc::SaveModified()
{
   CommitPendingEdit(); // so the "save changes?" question knows about it
   return __super::SaveModified();
}

BOOL CMappingEditorDoc::OnNewDocument()
{
   if (!CEAFDocument::OnNewDocument())
      return FALSE;

   m_Table = nlohmann::ordered_json::object();
   m_Table["format"] = "PGSuperIfcMapping";
   m_Table["version"] = 1;
   m_Table["name"] = "New table";
   m_Table["extends"] = "standard";
   LoadBaseTable();
   return TRUE;
}

BOOL CMappingEditorDoc::OpenTheDocument(LPCTSTR lpszPathName)
{
   std::ifstream file(std::filesystem::path(lpszPathName), std::ios::binary);
   if (!file)
   {
      CString strMsg;
      strMsg.Format(_T("%s can't be read."), lpszPathName);
      AfxMessageBox(strMsg, MB_OK | MB_ICONEXCLAMATION);
      return FALSE;
   }

   std::stringstream content;
   content << file.rdbuf();

   try
   {
      m_Table = nlohmann::ordered_json::parse(content.str());
   }
   catch (const std::exception& e)
   {
      CString strMsg;
      strMsg.Format(_T("%s isn't valid JSON, so the editor can't open it. Correct it in a text editor.\n\n%s"), lpszPathName, Utf8ToCString(e.what()).GetString());
      AfxMessageBox(strMsg, MB_OK | MB_ICONEXCLAMATION);
      return FALSE;
   }

   if (!m_Table.is_object() || m_Table.value("format", std::string()) != "PGSuperIfcMapping")
   {
      CString strMsg;
      strMsg.Format(_T("%s isn't an IFC mapping table (it doesn't have \"format\": \"PGSuperIfcMapping\")."), lpszPathName);
      AfxMessageBox(strMsg, MB_OK | MB_ICONEXCLAMATION);
      return FALSE;
   }

   SetPathName(lpszPathName, FALSE); // the base table is relative to the table's folder
   LoadBaseTable();
   return TRUE;
}

BOOL CMappingEditorDoc::SaveTheDocument(LPCTSTR lpszPathName)
{
   CommitPendingEdit();

   std::ofstream file(std::filesystem::path(lpszPathName), std::ios::binary);
   if (!file)
   {
      CString strMsg;
      strMsg.Format(_T("%s can't be written."), lpszPathName);
      AfxMessageBox(strMsg, MB_OK | MB_ICONEXCLAMATION);
      return FALSE;
   }

   file << FormatMappingTable(m_Table);
   return file.good() ? TRUE : FALSE;
}

void CMappingEditorDoc::DoIntegrateWithUI(BOOL bIntegrate)
{
   __super::DoIntegrateWithUI(bIntegrate);

   CEAFMainFrame* pFrame = EAFGetMainFrame();
   if (bIntegrate)
   {
      auto pStatusBar = new CMappingEditorStatusBar;
      pStatusBar->Create(pFrame);
      pFrame->SetStatusBar(pStatusBar); // the frame owns it
   }
   else
   {
      pFrame->SetStatusBar(nullptr); // back to the default status bar
   }
}

void CMappingEditorDoc::LoadDocumentSettings()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   __super::LoadDocumentSettings();
}

void CMappingEditorDoc::SaveDocumentSettings()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   __super::SaveDocumentSettings();
}

BOOL CMappingEditorDoc::GetStatusBarMessageString(UINT nID, CString& rMessage) const
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   return __super::GetStatusBarMessageString(nID, rMessage);
}

BOOL CMappingEditorDoc::GetToolTipMessageString(UINT nID, CString& rMessage) const
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   return __super::GetToolTipMessageString(nID, rMessage);
}

CString CMappingEditorDoc::GetToolbarSectionName()
{
   return _T("IfcMappingEditor");
}

CString CMappingEditorDoc::GetDocumentationRootLocation()
{
   return EAFGetApp()->GetDocumentationRootLocation();
}

HINSTANCE CMappingEditorDoc::GetResourceInstance()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   return AfxGetInstanceHandle();
}

void CMappingEditorDoc::DeleteContents()
{
   m_Table = nlohmann::ordered_json::object();
   m_pBase.reset();
   m_BaseError.clear();
   __super::DeleteContents();
}

void CMappingEditorDoc::OnValidate()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());

   CommitPendingEdit();

   std::string message;
   bool bValid = Validate(message);
   CMappingMessagesDlg dlg(bValid ? _T("Table Is Valid") : _T("Table Has Problems"), Utf8ToCString(message), EAFGetMainFrame());
   dlg.DoModal();
}
