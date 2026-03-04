///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright © 1999-2025  Washington State Department of Transportation
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

#include <Eigen/Dense>

struct Wire {
   std::vector<Eigen::Vector3d> verts;          // vertex positions
   std::vector<std::array<int, 2>> edges;        // edge list - edge[i] is connected by verts[edge[i][0]] and verts[edge[i][1]]
   std::vector<std::vector<int>> edge_edge;            // edge-edge adjacency - ee[i]

   Wire(const std::vector<Eigen::Vector3d>& v,
      const std::vector<std::array<int, 2>>& e);

   // ------------------------------------------------------------
   // Build adjacency: ee[i] = list of edges sharing a vertex with edge i
   // ------------------------------------------------------------
   void build_edge_adjacency();

   // ------------------------------------------------------------
   // Decompose wire into smooth-connected components
   // ------------------------------------------------------------
   std::vector<Wire> decompose(double max_edge_angle_deg = 10.0) const;

   // ------------------------------------------------------------
   // Total length of the wire - sum the length of each edge in the wire
   // ------------------------------------------------------------
   double length() const;

   // ------------------------------------------------------------
   // Bounding box center
   // ------------------------------------------------------------
   Eigen::Vector3d bbox_center() const;
};

struct Mesh {
   std::vector<Eigen::Vector3d> verts;              // vertex positions
   std::vector<std::array<int, 3>> faces;            // triangle indices
   std::map<std::pair<int, int>, std::vector<int>> edge2faces; // faces sharing an edge
   std::vector<std::vector<int>> face_face;         // face adjacency
   std::vector<Eigen::Vector3d> face_normals;       // per-face normals

   Mesh(const std::vector<Eigen::Vector3d>& v, const std::vector<std::array<int, 3>>& f);

   // ------------------------------------------------------------
   // Build edge to face mapping: edge2face[i] = list of faces sharing an edge
   // ------------------------------------------------------------
   void build_edge_map();


   // ------------------------------------------------------------
   // Build adjacency: face_face[i] = list of faces sharing an edge
   // ------------------------------------------------------------
   void build_adjacency();

   // ------------------------------------------------------------
   // Compute per-face normals (normalized)
   // ------------------------------------------------------------
   void compute_face_normals();

   // ------------------------------------------------------------
   // Decompose into smooth-connected components
   // ------------------------------------------------------------
   std::vector<Mesh> decompose(double max_edge_angle_deg = 10.0) const;

   // ------------------------------------------------------------
   // Boundary extraction (outer boundary only)
   // ------------------------------------------------------------
   Wire boundary() const;
};

Mesh get_top_mesh(const std::vector<Mesh>& meshes);