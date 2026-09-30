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
#include "IfcTableIds.h"

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
// CMappingFromIdsDlg

BEGIN_MESSAGE_MAP(CMappingFromIdsDlg, CDialog)
   ON_BN_CLICKED(IDC_FROMIDS_IDS_BROWSE, &CMappingFromIdsDlg::OnBrowseIds)
   ON_BN_CLICKED(IDC_FROMIDS_BINDING_BROWSE, &CMappingFromIdsDlg::OnBrowseBinding)
END_MESSAGE_MAP()

void CMappingFromIdsDlg::DoDataExchange(CDataExchange* pDX)
{
   CDialog::DoDataExchange(pDX);
   DDX_Text(pDX, IDC_FROMIDS_IDS, m_strIds);
   DDX_Text(pDX, IDC_FROMIDS_BINDING, m_strBinding);
   DDX_Check(pDX, IDC_FROMIDS_EXTEND, m_bExtend);
}

void CMappingFromIdsDlg::OnOK()
{
   if (!UpdateData(TRUE))
      return;
   m_strIds.Trim();
   m_strBinding.Trim();
   if (m_strIds.IsEmpty())
   {
      AfxMessageBox(_T("Choose the IDS."), MB_OK | MB_ICONEXCLAMATION);
      return;
   }
   CDialog::OnOK();
}

void CMappingFromIdsDlg::OnBrowseIds()
{
   CFileDialog dlg(TRUE, _T("ids"), nullptr, OFN_HIDEREADONLY | OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST, _T("IDS Files (*.ids)|*.ids|All Files (*.*)|*.*||"), this);
   if (dlg.DoModal() == IDOK)
      SetDlgItemText(IDC_FROMIDS_IDS, dlg.GetPathName());
}

void CMappingFromIdsDlg::OnBrowseBinding()
{
   CFileDialog dlg(TRUE, _T("json"), nullptr, OFN_HIDEREADONLY | OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST, _T("Binding Files (*.json)|*.json|All Files (*.*)|*.*||"), this);
   if (dlg.DoModal() == IDOK)
      SetDlgItemText(IDC_FROMIDS_BINDING, dlg.GetPathName());
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
   ON_COMMAND(ID_MAPPING_WRITE_IDS, &CMappingEditorDoc::OnWriteIds)
   ON_COMMAND(ID_MAPPING_FROM_IDS, &CMappingEditorDoc::OnGenerateFromIds)
   ON_COMMAND(ID_MAPPING_OPEN_MODEL, &CMappingEditorDoc::OnOpenModel)
   ON_COMMAND(ID_MAPPING_TRY_MODEL, &CMappingEditorDoc::OnTryModel)
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

const CIfcMappingTable* CMappingEditorDoc::GetCurrentTable(std::string& error)
{
   std::string text = FormatMappingTable(m_Table);
   if (!m_pCurrent || text != m_CurrentText)
   {
      m_CurrentText = text;
      m_CurrentError.clear();
      m_pCurrent.reset();
      try
      {
         m_pCurrent = CIfcMappingTable::Load(GetTablePath(), MappingTableSource::Editor, &text);
      }
      catch (const std::exception& e)
      {
         m_CurrentError = e.what();
      }
   }
   error = m_CurrentError;
   return m_pCurrent.get();
}

bool CMappingEditorDoc::OpenModel()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   CFileDialog dlg(TRUE, _T("ifc"), nullptr, OFN_HIDEREADONLY | OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST,
      _T("IFC Models (*.ifc)|*.ifc|All Files (*.*)|*.*||"), EAFGetMainFrame());
   if (dlg.DoModal() != IDOK)
      return false;

   try
   {
      CWaitCursor wait;
      m_pModel = CIfcModel::Open(std::filesystem::path(dlg.GetPathName().GetString()));
   }
   catch (const std::exception& e)
   {
      AfxMessageBox(Utf8ToCString(e.what()), MB_OK | MB_ICONEXCLAMATION);
      return false;
   }

   UpdateAllViews(nullptr, HINT_STRUCTURE); // the Model item shows the file
   SelectNode(MappingNode{ MappingNode::Kind::Model, "" });
   return true;
}

void CMappingEditorDoc::OnOpenModel()
{
   CommitPendingEdit();
   OpenModel();
}

void CMappingEditorDoc::OnTryModel()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   CommitPendingEdit();
   if (!m_pModel && !OpenModel())
      return;

   std::string error;
   const auto* pTable = GetCurrentTable(error);
   CString strReport;
   if (!pTable)
   {
      strReport = _T("The table must be valid to try it on a model.\n\n") + Utf8ToCString(error);
   }
   else
   {
      CWaitCursor wait;
      strReport = Utf8ToCString(m_pModel->TryTable(*pTable));
   }
   CMappingMessagesDlg dlg(_T("Try the Table on the Model"), strReport, EAFGetMainFrame());
   dlg.DoModal();
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
   m_pCurrent.reset();
   m_CurrentText.clear();
   m_CurrentError.clear();
   // the model stays open: it's the editor's, not the table's (e.g. when a table is generated from an IDS)
   __super::DeleteContents();
}

void CMappingEditorDoc::OnWriteIds()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   CommitPendingEdit();

   // the table as imports and exports see it
   std::unique_ptr<CIfcMappingTable> pTable;
   std::string text = FormatMappingTable(m_Table);
   try
   {
      pTable = CIfcMappingTable::Load(GetTablePath(), MappingTableSource::Editor, &text);
   }
   catch (const std::exception& e)
   {
      CMappingMessagesDlg dlg(_T("Table Has Problems"), _T("The table must be valid to write it as an IDS.\n\n") + Utf8ToCString(e.what()), EAFGetMainFrame());
      dlg.DoModal();
      return;
   }

   CString strDefault = GetPathName().IsEmpty() ? CString(_T("Untitled.ids")) : CString(GetTablePath().replace_extension(L".ids").filename().c_str());
   CFileDialog dlg(FALSE, _T("ids"), strDefault, OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, _T("IDS Files (*.ids)|*.ids|All Files (*.*)|*.*||"), EAFGetMainFrame());
   if (dlg.DoModal() != IDOK)
      return;

   std::ostringstream os;
   std::vector<std::string> notes;
   try
   {
      CWaitCursor wait;
      std::ofstream ids(std::filesystem::path(dlg.GetPathName().GetString()), std::ios::binary);
      if (!ids)
         throw std::runtime_error("The IDS file can't be written.");
      WriteTableAsIds(*pTable, ids, notes);
      os << "IDS written: " << CStringToUtf8(dlg.GetPathName()) << std::endl;
      for (const auto& note : notes)
         os << std::endl << "Note: " << note;
   }
   catch (const std::exception& e)
   {
      os << "The IDS wasn't written: " << e.what();
   }
   CMappingMessagesDlg result(_T("General IDS"), Utf8ToCString(os.str()), EAFGetMainFrame());
   result.DoModal();
}

void CMappingEditorDoc::OnGenerateFromIds()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   CommitPendingEdit();
   if (!SaveModified()) // the generated table replaces this one
      return;

   CMappingFromIdsDlg dlg(EAFGetMainFrame());
   if (dlg.DoModal() != IDOK)
      return;

   std::ostringstream os;
   try
   {
      CWaitCursor wait;
      // targets come from the standard table where the IDS and the binding file don't say
      auto standard = CIfcMappingTable::Load(std::filesystem::path(), MappingTableSource::InstalledStandard);
      auto result = GenerateTableFromIds(std::filesystem::path(dlg.m_strIds.GetString()), std::filesystem::path(dlg.m_strBinding.GetString()), *standard, "", dlg.m_bExtend ? true : false);

      m_Table = nlohmann::ordered_json::parse(result.table_json);

      // a new table, to be saved
      m_strPathName.Empty();
      SetTitle(_T("Untitled"));
      for (POSITION pos = GetFirstViewPosition(); pos; )
      {
         if (CFrameWnd* pFrame = GetNextView(pos)->GetParentFrame())
            pFrame->OnUpdateFrameTitle(TRUE);
      }
      static_cast<CFrameWnd*>(EAFGetMainFrame())->OnUpdateFrameTitle(TRUE); // public in CFrameWnd, virtual
      LoadBaseTable();
      SetModifiedFlag(TRUE);
      UpdateAllViews(nullptr, HINT_STRUCTURE);
      SelectNode(MappingNode{ MappingNode::Kind::Table, "" });

      os << "Table generated from " << CStringToUtf8(dlg.m_strIds);
      if (!dlg.m_strBinding.IsEmpty())
         os << " with the binding file " << CStringToUtf8(dlg.m_strBinding);
      os << "." << std::endl << result.bound << " property facets with a target, " << result.unbound << " without." << std::endl;
      for (const auto& line : result.report)
         os << std::endl << line;

      std::string message;
      bool bValid = Validate(message);
      os << std::endl << std::endl << (bValid ? "The generated table is valid." : "The generated table has problems:\n" + message);
   }
   catch (const std::exception& e)
   {
      os << "The table wasn't generated: " << e.what();
   }
   CMappingMessagesDlg result(_T("Table Generated from IDS"), Utf8ToCString(os.str()), EAFGetMainFrame());
   result.DoModal();
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
