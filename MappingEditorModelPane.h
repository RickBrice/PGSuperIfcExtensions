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

// Mapping table editor, stage 4: the Model item. Browse the elements of an element role in an IFC model and their
// properties, see what the table reads for a target and the importer's hints, and add a property as a location.

#include "MappingEditorPanes.h"
#include "IfcModelBrowser.h"

class CModelPane : public CMappingPane
{
public:
   CModelPane(CMappingEditorDoc* pDoc) : CMappingPane(IDD_MAPPING_MODEL_PANE, pDoc) {}

protected:
   BOOL OnInitDialog() override;
   afx_msg void OnOpenModel();
   afx_msg void OnRoleChanged();
   afx_msg void OnElementChanged();
   afx_msg void OnTargetChanged();
   afx_msg void OnAdd();
   afx_msg void OnPropertiesChanged(NMHDR* pNMHDR, LRESULT* pResult);
   afx_msg void OnPropertiesDblClk(NMHDR* pNMHDR, LRESULT* pResult);
   DECLARE_MESSAGE_MAP()

private:
   CListCtrl m_Properties;
   std::vector<ElementKind> m_Roles;
   std::vector<ModelElement> m_Elements;
   std::vector<ModelProperty> m_ElementProperties;
   std::vector<std::string> m_Targets;

   const CIfcMappingTable* GetTable(std::string* pError = nullptr);
   int GetElementId();
   std::string GetTargetName();
   void UpdateReading();
   void UpdateButtons();
};
