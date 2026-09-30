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
#include <EAF\EAFStatusBar.h>
#include <EAF\EAFDocTemplate.h>
#include <EAF\ComponentObject.h>
#include <EAF\PluginApp.h>
#include <nlohmann/json.hpp>

#include "IfcMappingTable.h"

#include <filesystem>
#include <memory>
#include <string>

// UTF-8 (the table's text) <-> the UI's text
CString Utf8ToCString(const std::string& text);
std::string CStringToUtf8(const CString& text);

// The editor's status bar. It must have the AutoSave indicator: the data recovery handler writes to it
// (CEAFStatusBar::AutoSaveSaving), and the default status bar doesn't have one
class CMappingEditorStatusBar : public CEAFStatusBar
{
protected:
   void GetStatusIndicators(const UINT** lppIDArray, int* pnIDCount) override;
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

class CMappingEditorDoc : public CEAFDocument
{
protected:
   CMappingEditorDoc();
   DECLARE_DYNCREATE(CMappingEditorDoc)

public:
   // UpdateAllViews hints
   static constexpr LPARAM HINT_STRUCTURE = 1; // what the table defines changed (tree labels)
   static constexpr LPARAM HINT_CONTENT = 2;   // values changed

   nlohmann::ordered_json& GetTable() { return m_Table; }

   // Called by the panes after they change the table
   void TableChanged(bool bStructure);

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

   BOOL OnNewDocument() override;
   BOOL SaveModified() override;
   BOOL OpenTheDocument(LPCTSTR lpszPathName) override;
   BOOL SaveTheDocument(LPCTSTR lpszPathName) override;

   void DoIntegrateWithUI(BOOL bIntegrate) override;
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
   DECLARE_MESSAGE_MAP()

private:
   nlohmann::ordered_json m_Table;
   std::unique_ptr<CIfcMappingTable> m_pBase;
   std::string m_BaseError;
};
