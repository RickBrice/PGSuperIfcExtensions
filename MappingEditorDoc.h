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

/*****************************************************************************
   Mapping table editor (devdocs/IfcImportPlan.md, R3; M7)

   A BridgeLink application in the IFC extension DLL that edits one mapping
   table file at a time. The document is the table's JSON (an ordered
   nlohmann JSON tree), so key order, comments, and anything the editor
   doesn't show are kept. Validation loads the JSON with CIfcMappingTable,
   as imports and exports do, so a table the editor accepts is a table the
   importer accepts. Files are written in the style of Standard.json
   (FormatMappingTable).
*****************************************************************************/

#include <EAF\EAFDocument.h>
#include <EAF\EAFDocTemplate.h>
#include <EAF\ComponentObject.h>
#include <EAF\PluginApp.h>
#include <nlohmann/json.hpp>

#include "IfcMappingTable.h"
#include "IfcModelBrowser.h"

#include <filesystem>
#include <memory>
#include <string>

// UTF-8 (the table's text) <-> the UI's text
CString Utf8ToCString(const std::string& text);
std::string CStringToUtf8(const CString& text);

// Table > Generate from IDS: the IDS, the binding file, and whether the table extends the standard table
class CMappingFromIdsDlg : public CDialog
{
public:
   CMappingFromIdsDlg(CWnd* pParent = nullptr) : CDialog(IDD_MAPPING_FROM_IDS, pParent) {}

   CString m_strIds;
   CString m_strBinding;
   BOOL m_bExtend = TRUE;

protected:
   void DoDataExchange(CDataExchange* pDX) override;
   void OnOK() override;
   afx_msg void OnBrowseIds();
   afx_msg void OnBrowseBinding();
   DECLARE_MESSAGE_MAP()
};

// The BridgeLink application plugin
class CMappingEditorPlugin : public WBFL::EAF::ComponentObject, public WBFL::EAF::IPluginApp
{
public:
   BOOL Init(CEAFApp* pParent) override;
   void Terminate() override;
   void IntegrateWithUI(BOOL bIntegrate) override;
   std::vector<CEAFDocTemplate*> CreateDocTemplates() override;
   HMENU GetSharedMenuHandle() override;
   CString GetName() override;
   CString GetDocumentationSetName() override;
   CString GetDocumentationURL() override;
   CString GetDocumentationMapFile() override;
   void LoadDocumentationMap() override;
   std::pair<WBFL::EAF::HelpResult, CString> GetDocumentLocation(LPCTSTR lpszDocSetName, UINT nID) override;
};

class CMappingEditorDocTemplate : public CEAFDocTemplate
{
public:
   CMappingEditorDocTemplate(UINT nIDResource, std::shared_ptr<WBFL::EAF::ICommandCallback> pCallback, CRuntimeClass* pDocClass,
      CRuntimeClass* pFrameClass, CRuntimeClass* pViewClass, HMENU hSharedMenu = nullptr, int maxViewCount = -1);

   CString GetTemplateGroupItemDescription(const CEAFTemplateItem* pItem) const override;

   DECLARE_DYNAMIC(CMappingEditorDocTemplate)
};

struct MappingNode;

class CMappingEditorDoc : public CEAFDocument
{
protected:
   CMappingEditorDoc();
   DECLARE_DYNCREATE(CMappingEditorDoc)

public:
   // UpdateAllViews hints
   static constexpr LPARAM HINT_STRUCTURE = 1; // what the table defines changed (tree labels)
   static constexpr LPARAM HINT_CONTENT = 2;   // values changed
   static constexpr LPARAM HINT_SELECT = 3;    // select a tree node (CMappingSelectHint), after the tree is rebuilt

   nlohmann::ordered_json& GetTable() { return m_Table; }

   // Called by the panes after they change the table
   void TableChanged(bool bStructure);

   // Selects a tree node and shows its pane, once the current message is handled (so a pane can ask for it and be replaced)
   void SelectNode(const MappingNode& node);

   // The tables named by "extends" (nullptr if none, or if they can't be loaded; see GetBaseTableError)
   const CIfcMappingTable* GetBaseTable() const { return m_pBase.get(); }
   const std::string& GetBaseTableError() const { return m_BaseError; }
   void LoadBaseTable();

   // The table's file, or Untitled.json in the current folder for a new table (a relative "extends" is relative to its folder)
   std::filesystem::path GetTablePath() const;

   // Loads the table as imports and exports do. message: the table and the tables it extends, or the error
   bool Validate(std::string& message) const;

   // Writes a field being edited to the table (a pane writes a field when it loses the focus)
   void CommitPendingEdit();

   // The table as imports and exports see it: loaded from the current JSON, again when the JSON changes.
   // nullptr, with the table's error, if it isn't valid
   const CIfcMappingTable* GetCurrentTable(std::string& error);

   // The IFC model open in the editor (stage 4), or nullptr
   CIfcModel* GetModel() { return m_pModel.get(); }
   const CIfcModel* GetModel() const { return m_pModel.get(); }

   // Asks for a model file and opens it, and shows the Model item. Returns true if a model was opened
   bool OpenModel();

   BOOL OnNewDocument() override;
   BOOL SaveModified() override;
   BOOL OpenTheDocument(LPCTSTR lpszPathName) override;
   BOOL SaveTheDocument(LPCTSTR lpszPathName) override;

   void LoadDocumentSettings() override;
   void SaveDocumentSettings() override;
   BOOL GetStatusBarMessageString(UINT nID, CString& rMessage) const override;
   BOOL GetToolTipMessageString(UINT nID, CString& rMessage) const override;
   CString GetToolbarSectionName() override;
   CString GetDocumentationRootLocation() override;

protected:
   HINSTANCE GetResourceInstance() override;
   void DeleteContents() override;

   afx_msg void OnValidate();
   afx_msg void OnWriteIds();
   afx_msg void OnGenerateFromIds();
   afx_msg void OnOpenModel();
   afx_msg void OnTryModel();
   DECLARE_MESSAGE_MAP()

private:
   nlohmann::ordered_json m_Table;
   std::unique_ptr<CIfcMappingTable> m_pBase;
   std::string m_BaseError;

   std::unique_ptr<CIfcModel> m_pModel;
   std::unique_ptr<CIfcMappingTable> m_pCurrent; // GetCurrentTable
   std::string m_CurrentText;                    // the JSON m_pCurrent was loaded from
   std::string m_CurrentError;
};
