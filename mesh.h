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

  Vector3T<Real> nml; // face normal
private:
  Size v[3];
};


class Vertex
{
public:
  explicit Vertex() : pos(0.0f,0.0f,0.0f), nml(0.0f,0.0f,0.0f) { }
  Vector3T<Real> pos; // vertex position
  Vector3T<Real> nml; // vertex normal
};

typedef std::vector<Vertex> VertexVec;
typedef std::vector<Face> FaceVec;

// A triangular mesh.

class Mesh : public Primitive 
{
public:
  explicit Mesh(const aiMesh *mesh, Index matidx);
  explicit Mesh(const Mesh &orig);
  virtual ~Mesh();

  AlignedBox3f &compute_bbox();
  void compute_face_normals();

  void transform_in_place(const Affine3f &trans);

  inline const VertexVec &get_verts() const { return m_verts; }
  inline const FaceVec   &get_faces() const { return m_faces; }
  inline bool is_mesh() const { return true; }

  // manage position and normal visualization data
  void prepare_visposnml();
  inline const Matrix3Xf &get_vispos() const { return m_vispos; }
  inline const Matrix3Xf &get_visnml() const { return m_visnml; }
  inline bool is_staleposnml() { return m_staleposnml; }
  inline void set_staleposnml(bool sp) { m_staleposnml = sp; }

  friend std::ostream& operator<<(std::ostream& out, const Mesh& mesh);

protected:
  VertexVec m_verts;
  FaceVec   m_faces;

  // The following data structures provide a mechanism for syncronizing computed
  // mesh data with its visual representative (like glmesh)
  Matrix3Xf m_vispos;
  Matrix3Xf m_visnml;
  std::atomic<bool> m_staleposnml;
};

typedef std::vector< Mesh >       MeshVec;
typedef boost::shared_ptr< Mesh > MeshPtr;

#endif // MESH_H
