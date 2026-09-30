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
#include "stdafx.h"
#include "IfcTableFormat.h"

#include <sstream>
#include <vector>

using ordered_json = nlohmann::ordered_json;

namespace
{
   const char* const KEY_ORDER[] = { "comment", "format", "version", "name", "extends", "elements", "property_sets", "quantity_sets",
      "classification_systems", "classifications", "targets" };

   // A value on one line, with ", " and ": " separators and UTF-8 text as it is (as Python's json.dumps(ensure_ascii=False))
   std::string inline_json(const ordered_json& v)
   {
      if (v.is_object())
      {
         std::string s = "{";
         bool first = true;
         for (const auto& [key, value] : v.items())
         {
            s += (first ? "" : ", ") + ordered_json(key).dump() + ": " + inline_json(value);
            first = false;
         }
         return s + "}";
      }
      if (v.is_array())
      {
         std::string s = "[";
         for (size_t i = 0; i < v.size(); i++)
            s += (i == 0 ? "" : ", ") + inline_json(v[i]);
         return s + "]";
      }
      return v.dump();
   }

   std::string quote(const std::string& key)
   {
      return ordered_json(key).dump();
   }

   // one item per line
   void write_list(std::vector<std::string>& lines, const ordered_json& list, const std::string& indent)
   {
      for (size_t i = 0; i < list.size(); i++)
         lines.push_back(indent + inline_json(list[i]) + (i + 1 < list.size() ? "," : ""));
   }

   void write_section(std::vector<std::string>& lines, const std::string& key, const ordered_json& v, const std::string& end)
   {
      if (key == "elements" && v.is_object())
      {
         lines.push_back("  " + quote(key) + ": {");
         size_t n = 0;
         for (const auto& [role, selector] : v.items())
            lines.push_back("    " + quote(role) + ": " + inline_json(selector) + (++n < v.size() ? "," : ""));
         lines.push_back("  }" + end);
      }
      else if (key == "targets" && v.is_object())
      {
         lines.push_back("  " + quote(key) + ": {");
         size_t n = 0;
         for (const auto& [target, item] : v.items())
         {
            std::string sep = (++n < v.size()) ? "," : "";
            if (item.is_array())
            {
               lines.push_back("    " + quote(target) + ": [");
               write_list(lines, item, "      ");
               lines.push_back("    ]" + sep);
            }
            else if (item.is_object() && item.contains("locations") && item["locations"].is_array())
            {
               // { "mode": "replace", "locations": [...] }
               std::string head = "    " + quote(target) + ": {";
               for (const auto& [k, value] : item.items())
               {
                  if (k != "locations")
                     head += quote(k) + ": " + inline_json(value) + ", ";
               }
               lines.push_back(head + "\"locations\": [");
               write_list(lines, item["locations"], "      ");
               lines.push_back("    ]}" + sep);
            }
            else
            {
               lines.push_back("    " + quote(target) + ": " + inline_json(item) + sep);
            }
         }
         lines.push_back("  }" + end);
      }
      else if ((key == "property_sets" || key == "quantity_sets") && v.is_array())
      {
         const std::string items_key = (key == "property_sets") ? "properties" : "quantities";
         lines.push_back("  " + quote(key) + ": [");
         for (size_t n = 0; n < v.size(); n++)
         {
            const auto& set = v[n];
            std::string sep = (n + 1 < v.size()) ? "," : "";
            if (!set.is_object() || !set.contains(items_key) || !set[items_key].is_array())
            {
               lines.push_back("    " + inline_json(set) + sep);
               continue;
            }

            lines.push_back("    {");
            std::string head;
            for (const auto& [k, value] : set.items())
            {
               if (k != items_key && k != "comment")
                  head += (head.empty() ? "" : ", ") + quote(k) + ": " + inline_json(value);
            }
            if (!head.empty())
               lines.push_back("      " + head + ",");
            if (set.contains("comment"))
               lines.push_back("      \"comment\": " + inline_json(set["comment"]) + ",");
            lines.push_back("      " + quote(items_key) + ": [");
            write_list(lines, set[items_key], "        ");
            lines.push_back("      ]");
            lines.push_back("    }" + sep);
         }
         lines.push_back("  ]" + end);
      }
      else if ((key == "classification_systems" || key == "classifications") && v.is_array())
      {
         lines.push_back("  " + quote(key) + ": [");
         write_list(lines, v, "    ");
         lines.push_back("  ]" + end);
      }
      else
      {
         lines.push_back("  " + quote(key) + ": " + inline_json(v) + end);
      }
   }
}

std::string FormatMappingTable(const ordered_json& table)
{
   // the sections in their order, then any other keys in the order they appear
   std::vector<std::string> keys;
   for (const char* key : KEY_ORDER)
   {
      if (table.contains(key))
         keys.push_back(key);
   }
   for (const auto& [key, value] : table.items())
   {
      if (std::find(keys.begin(), keys.end(), key) == keys.end())
         keys.push_back(key);
   }

   std::vector<std::string> lines{ "{" };
   for (size_t i = 0; i < keys.size(); i++)
      write_section(lines, keys[i], table[keys[i]], i + 1 < keys.size() ? "," : "");
   lines.push_back("}");

   std::string text;
   for (const auto& line : lines)
      text += line + "\r\n";
   return text;
}
