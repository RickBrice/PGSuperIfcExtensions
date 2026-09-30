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

// Helpers shared by the mapping table editor's panes and dialogs

#include "MappingEditorDoc.h"
#include "IfcTargets.h"
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

namespace mapping_editor
{
   using ordered_json = nlohmann::ordered_json;

   inline CString GetText(CWnd* pWnd, int nID)
   {
      CString text;
      pWnd->GetDlgItemText(nID, text);
      return text;
   }

   // Multi-line edit controls need CR LF
   inline CString ToWindowsLines(CString text)
   {
      text.Replace(_T("\r\n"), _T("\n"));
      text.Replace(_T("\n"), _T("\r\n"));
      return text;
   }

   // The lines of a multi-line edit control, trimmed, without empty lines
   inline std::vector<CString> GetLines(CWnd* pWnd, int nID)
   {
      std::vector<CString> lines;
      CString text = GetText(pWnd, nID);
      int pos = 0;
      for (CString line = text.Tokenize(_T("\r\n"), pos); !line.IsEmpty(); line = text.Tokenize(_T("\r\n"), pos))
      {
         line.Trim();
         if (!line.IsEmpty())
            lines.push_back(line);
      }
      return lines;
   }

   inline std::string string_value(const ordered_json& j, const char* key)
   {
      return (j.is_object() && j.contains(key) && j[key].is_string()) ? j[key].get<std::string>() : std::string();
   }

   // Sets a text member, or removes it when the text is empty. Returns true if the object changed
   inline bool set_or_erase(ordered_json& j, const char* key, const std::string& value)
   {
      if (value.empty())
      {
         if (!j.contains(key))
            return false;
         j.erase(key);
         return true;
      }
      if (j.contains(key) && j[key] == value)
         return false;
      j[key] = value;
      return true;
   }

   // A one-line JSON value, as the table file shows it
   inline std::string inline_json(const ordered_json& j)
   {
      return j.dump(-1, ' ', false);
   }

   // The element roles of "applies_to" (a role or a list of roles)
   inline std::vector<std::string> roles_of(const ordered_json& item)
   {
      std::vector<std::string> roles;
      if (!item.is_object() || !item.contains("applies_to"))
         return roles;
      const auto& applies_to = item["applies_to"];
      if (applies_to.is_string())
         roles.push_back(applies_to.get<std::string>());
      else if (applies_to.is_array())
      {
         for (const auto& role : applies_to)
         {
            if (role.is_string())
               roles.push_back(role.get<std::string>());
         }
      }
      return roles;
   }

   // All element role names, in ElementKind order
   inline std::vector<std::string> all_roles()
   {
      std::vector<std::string> roles;
      for (int kind = (int)ElementKind::Project; kind <= (int)ElementKind::Barrier; kind++)
         roles.emplace_back(GetElementRoleName((ElementKind)kind));
      return roles;
   }

   // Fills a multiple-selection list box with the element roles and selects the given ones
   inline void FillRoleList(CListBox* pList, const std::vector<std::string>& selected)
   {
      int first = -1;
      for (const auto& role : all_roles())
      {
         int idx = pList->AddString(Utf8ToCString(role));
         if (std::find(selected.begin(), selected.end(), role) != selected.end())
         {
            pList->SetSel(idx, TRUE);
            if (first < 0)
               first = idx;
         }
      }
      if (0 <= first)
         pList->SetTopIndex(first); // so the selected role shows
   }

   // The selected roles of the list box as "applies_to": a role, or a list of roles
   inline ordered_json GetRoleList(CListBox* pList)
   {
      ordered_json roles = ordered_json::array();
      for (int i = 0; i < pList->GetCount(); i++)
      {
         if (0 < pList->GetSel(i))
         {
            CString text;
            pList->GetText(i, text);
            roles.push_back(CStringToUtf8(text));
         }
      }
      return roles.size() == 1 ? roles[0] : roles;
   }
}
