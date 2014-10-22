#ifndef MESH_H
#define MESH_H

#include <vector>
#include <iostream>
#include <assimp/mesh.h>
#include <boost/shared_ptr.hpp>
#include "types.h"
#include "primitive.h"


class Mesh;

std::ostream& operator<<(std::ostream& out, const Mesh& mesh);

class Face {
public:
  explicit Face() { }
  explicit Face(Size v0, Size v1, Size v2) { v[0] = v0; v[1] = v1; v[2] = v2; }
  Size  &operator[](Size i) { return v[i]; }
  Size   operator[](Size i) const { return v[i]; }

  Vector3R<Real> nml; // face normal
private:
  Size v[3];
};


class Vertex
{
public:
  explicit Vertex() : pos(0.0f,0.0f,0.0f), nml(0.0f,0.0f,0.0f) { }
  Vector3R<Real> pos; // vertex position
  Vector3R<Real> nml; // vertex normal
};

typedef std::vector<Vertex> VertexVec;
typedef std::vector<Face> FaceVec;

// A triangular mesh.

class Mesh : public Primitive 
{
public:
  explicit Mesh(const aiMesh *mesh);
  ~Mesh();

  AlignedBox3f &compute_bbox();
  void compute_face_normals();

  void transform_in_place(const AffineCompact3f &trans);

  inline const VertexVec &get_verts() const { return m_verts; }
  inline const FaceVec   &get_faces() const { return m_faces; }
  inline bool is_mesh() const { return true; }

  friend std::ostream& operator<<(std::ostream& out, const Mesh& mesh);

protected:
  VertexVec m_verts;
  FaceVec   m_faces;
};

typedef boost::shared_ptr< Mesh > MeshPtr;

#endif // MESH_H
