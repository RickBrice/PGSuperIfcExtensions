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

#include "Geometry.h"
#include "Utilities.h"

#include <stack>

Wire::Wire(const std::vector<Eigen::Vector3d>& v,
   const std::vector<std::array<int, 2>>& e)
   : verts(v), edges(e)
{
   build_edge_adjacency();
}

// ------------------------------------------------------------
// Build adjacency: ee[i] = list of edges sharing a vertex with edge i
// ------------------------------------------------------------
void Wire::build_edge_adjacency() {
   std::map<int, std::vector<int>> v2e; // key = vertex index. value = list of edges that connect to that vertex

   // Build vertex → edges map
   for (int ei = 0; ei < edges.size(); ++ei) {
      int a = edges[ei][0];
      int b = edges[ei][1];
      v2e[a].push_back(ei);
      v2e[b].push_back(ei);
   }

   // Build adjacency sets
   std::vector<std::set<int>> neigh(edges.size());
   for (auto& kv : v2e) {
      const auto& es = kv.second; // edges sharing a vertex
      if (es.size() > 1) {
         for (int i : es)
         {
            for (int j : es)
            {
               if (i != j)
                  neigh[i].insert(j); // edge j is a neighbor of edge i
            }
         }
      }
   }

   // convert to list of adjacent edge indices
   edge_edge.resize(edges.size());
   for (int i = 0; i < neigh.size(); ++i)
      edge_edge[i] = std::vector<int>(neigh[i].begin(), neigh[i].end());
}

// ------------------------------------------------------------
// Decompose wire into smooth-connected components
// ------------------------------------------------------------
std::vector<Wire> Wire::decompose(double max_edge_angle_deg) const {
   double cos_thr = std::cos(deg2rad(max_edge_angle_deg));

   std::vector<bool> used(edges.size(), false);
   std::vector<Wire> components;

   for (int seed = 0; seed < edges.size(); ++seed) {
      if (used[seed]) // if this edge already belongs to a sub-wire, skip it
         continue;

      // seed the stack with an edge, look for all other edges that are approximately in the same direction
      std::stack<int> stack;
      stack.push(seed);
      used[seed] = true;

      std::vector<int> comp; // indices of edges that are part of the current sub-wire

      while (!stack.empty()) {
         int i = stack.top();
         stack.pop();
         comp.push_back(i);

         for (int j : edge_edge[i]) {
            if (used[j])
               continue; // this edge is already used in a different sub-wire so skip it

            // Compute direction vectors
            Eigen::Vector3d v1 = verts[edges[i][1]] - verts[edges[i][0]];
            Eigen::Vector3d v2 = verts[edges[j][1]] - verts[edges[j][0]];

            double n1 = v1.norm();
            double n2 = v2.norm();
            if (n1 < 1e-12 || n2 < 1e-12)
               continue;

            v1 /= n1; // normalize vectors
            v2 /= n2;

            if (v1.dot(v2) >= cos_thr) {
               used[j] = true; // edges hare about in the same direction
               stack.push(j); // add edge index to the stack so the neighboring edges of this edge get checked next time through the while loop
            }
         }
      }

      // Build sub-wire
      std::vector<std::array<int, 2>> sub_edges;
      sub_edges.reserve(comp.size());
      for (int ei : comp)
         sub_edges.push_back(edges[ei]);

      components.emplace_back(verts, sub_edges); // save the sub-wire
   }

   return components; // return all sub-wires
}

// ------------------------------------------------------------
// Total length of the wire - sum the length of each edge in the wire
// ------------------------------------------------------------
double Wire::length() const {
   double L = 0.0;
   for (auto& e : edges) {
      Eigen::Vector3d d = verts[e[1]] - verts[e[0]];
      L += d.norm();
   }
   return L;
}

double Wire::plan_length() const {
   double L = 0.0;
   for (auto& e : edges) {
      Eigen::Vector2d d = verts[e[1]].block<2, 1>(0,0) - verts[e[0]].block<2, 1>(0, 0);
      L += d.norm();
   }
   return L;
}

// ------------------------------------------------------------
// Bounding box center
// ------------------------------------------------------------
Eigen::Vector3d Wire::bbox_center() const {
   if (edges.empty())
      return Eigen::Vector3d::Zero();

   // create bound box around all the edges of the wire
   Eigen::Vector3d vmin = verts[edges[0][0]];
   Eigen::Vector3d vmax = vmin;

   for (auto& e : edges) {
      for (int vi : e) {
         vmin = vmin.cwiseMin(verts[vi]);
         vmax = vmax.cwiseMax(verts[vi]);
      }
   }

   // center point is average value
   return 0.5 * (vmin + vmax);
}

Mesh::Mesh(const std::vector<Eigen::Vector3d>& v, const std::vector<std::array<int, 3>>& f,double face_normal_scale) :
   verts(v), faces(f), face_normal_scale(face_normal_scale)
{
   build_edge_map();
   build_adjacency();
   compute_face_normals();
}

// ------------------------------------------------------------
// Build edge to face mapping: edge2face[i] = list of faces sharing an edge
// ------------------------------------------------------------
void Mesh::build_edge_map()
{
   for (int fi = 0; fi < faces.size(); ++fi) {
      auto [a, b, c] = faces[fi]; // the three vertex indices that make up this face

      for (auto [u, v] : { std::pair{a,b}, std::pair{b,c}, std::pair{c,a} }) {
         if (u > v) std::swap(u, v); // order the edges so (7,3) and (3,7) are treated the same
         edge2faces[{u, v}].push_back(fi); // for each edge (u,v) store the face index
      }
   }
}


// ------------------------------------------------------------
// Build adjacency: face_face[i] = list of faces sharing an edge
// ------------------------------------------------------------
void Mesh::build_adjacency() {
   // for each edge shared by multiple faces, the faces are considered to be neighbors
   // if edge (u,v) is used by faces m,n,o, then m, n, o are neighbors
   std::vector<std::set<int>> neigh(faces.size()); // automatic sorting since we are using std::set

   for (auto& kv : edge2faces) {
      const auto& fs = kv.second; // faces sharing an edge
      if (fs.size() > 1) {
         for (int i : fs)
         {
            for (int j : fs)
            {
               if (i != j)
                  neigh[i].insert(j); // face j is a neighbor of face i
            }
         }
      }
   }

   // convert to list of adjacent face indices
   face_face.resize(faces.size());
   for (int i = 0; i < neigh.size(); ++i)
      face_face[i] = std::vector<int>(neigh[i].begin(), neigh[i].end());
}

// ------------------------------------------------------------
// Compute per-face normals (normalized)
// ------------------------------------------------------------
void Mesh::compute_face_normals() {
   face_normals.resize(faces.size());

   double volume = 0.0;

   for (int i = 0; i < faces.size(); ++i) {
      auto [a, b, c] = faces[i];
      const Eigen::Vector3d& v0 = verts[a];
      const Eigen::Vector3d& v1 = verts[b];
      const Eigen::Vector3d& v2 = verts[c];

      Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);
      double len = n.norm();
      if (len > 1e-12)
         n /= len;

      face_normals[i] = n;

      // volume of tetrahedron formed by the face vertices
      double vol = v0.dot(v1.cross(v2)); // could divide by 6 here, but why since all in the sum need to be divided by 6, do it after the loop
      volume += vol;
   }

   // divide 6 here to get the volume
   // though, this isn't strictly needed because we only care about the sign
   // and dividing by 6 doesn't change that
   volume /= 6.0;

   // now determine if the mesh is open or closed. it is closed if all edges are bound by exactly two faces
   bool closed = true;
   for (auto& edge : edge2faces)
   {
      if (edge.second.size() != 2)
      {
         closed = false;
         break;
      }
   }

   if (closed)
   {
      // volume > 0, normals point outward
      // volume < 0, normals point inward
      // volume approx = 0, mesh is degenerate or not watertight
      if (volume < 0.0)
      {
         face_normal_scale = -1.0;
      }
   }
}

// ------------------------------------------------------------
// Decompose into smooth-connected components
// ------------------------------------------------------------
std::vector<Mesh> Mesh::decompose(double max_edge_angle_deg) const {
   double cos_thr = std::cos(deg2rad(max_edge_angle_deg));

   std::vector<bool> used(faces.size(), false);
   std::vector<Mesh> components;

   for (int seed = 0; seed < faces.size(); ++seed) {
      if (used[seed]) // if this face already belongs to a sub-mesh, skip it
         continue;

      // seed the stack with a face, look for all other faces that have approximately the same surface normal
      std::stack<int> stack;
      stack.push(seed);
      used[seed] = true;

      std::vector<int> comp; // indices of faces that are part of the current sub-mesh

      while (!stack.empty()) {
         int i = stack.top();
         stack.pop();
         comp.push_back(i);

         const Eigen::Vector3d& ni = face_normals[i]; // surface normal of the curren face

         // check all neighboring faces to see if they have the same surface normal
         for (int j : face_face[i]) {
            if (used[j])
               continue; // this face is already used in a different sub-mesh so skip it

            if (ni.dot(face_normals[j]) >= cos_thr) {
               used[j] = true; // faces have about the same surface normal
               stack.push(j); // add face index to the stack so the neighboring faces of this face get checked next time through the while loop
            }
         }
      }

      // Build sub-mesh
      std::vector<std::array<int, 3>> sub_faces; // vector of faces, defined by their vertex indices, which all have about the same surface normal
      sub_faces.reserve(comp.size());
      for (int fi : comp)
         sub_faces.push_back(faces[fi]);

      components.emplace_back(verts, sub_faces, face_normal_scale); // save the sub-mesh
   }

   return components; // return all sub-meshes
}

// ------------------------------------------------------------
// Boundary extraction (outer boundary only)
// ------------------------------------------------------------
Wire Mesh::boundary() const {
   // Collect boundary edges (edges used by exactly 1 face)
   std::vector<std::array<int, 2>> boundary_edges; // bound_edges[i] -> start vertex index, end vertex index
   for (auto& kv : edge2faces) {
      if (kv.second.size() == 1) {
         boundary_edges.push_back({ kv.first.first, kv.first.second });
      }
   }

   // Build a graph to order boundary edges into a cycle
   // (Assumes a single outer boundary)
   std::map<int, std::vector<int>> graph;
   for (auto& e : boundary_edges) {
      graph[e[0]].push_back(e[1]);
      graph[e[1]].push_back(e[0]);
   }

   // Find cycle
   std::vector<int> cycle;
   int start = boundary_edges[0][0];
   int prev = -1;
   int cur = start;

   do {
      cycle.push_back(cur);
      const auto& nbrs = graph[cur];
      int next = (nbrs[0] == prev ? nbrs[1] : nbrs[0]);
      prev = cur;
      cur = next;
   } while (cur != start);

   // Convert cycle to edges
   std::vector<std::array<int, 2>> ordered_edges;
   for (int i = 0; i < cycle.size(); ++i) {
      int a = cycle[i];
      int b = cycle[(i + 1) % cycle.size()];
      ordered_edges.push_back({ a, b });
   }

   return Wire(verts, ordered_edges);
}

Mesh get_top_mesh(const std::vector<Mesh>& meshes)
{
   auto top_component_it = std::max_element(
      meshes.begin(), meshes.end(),
      [](const Mesh& a, const Mesh& b)
      {
         auto score = [](const Mesh& c)
            {
               // Compute average normal
               Eigen::Vector3d avg = Eigen::Vector3d::Zero();
               for (auto& n : c.face_normals)
                  avg += n;

               avg /= (double)c.face_normals.size();
               avg *= c.face_normal_scale; // apply adjustment factor for possibly inverted meshes

               // Compute max Z of (vertex + avg_normal), this essentially explodes the mesh so the top is obvious
               double best = -std::numeric_limits<double>::infinity();
               for (auto& f : c.faces) {
                  for (int vi : f) {
                     double z = (c.verts[vi] + avg).z();
                     best = std::max(best,z);
                  }
               }
               return best;
            };

         return score(a) < score(b);
      }
   );

   return *top_component_it;
}

void Mesh::print(std::ostream& os) const
{
   os << "verts" << std::endl;
   int idx = 0;
   Eigen::IOFormat fmt(8, 0, ", ", ", ", "", "", "", "");
   for (auto& v : verts)
   {
      os << idx << ", " << v.format(fmt) << std::endl;
      idx++;
   }

   idx = 0;
   os << "faces" << std::endl;
   for (auto& f : faces)
   {
      os << idx << ", " << f[0] << ", " << f[1] << ", " << f[2] << std::endl;
      idx++;
   }
}
