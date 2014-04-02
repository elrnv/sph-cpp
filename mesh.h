#ifndef MESH_H
#define MESH_H

#include <vector>
#include <iostream>
#include <assimp/mesh.h>
#include <boost/shared_ptr.hpp>
#include "primitive.h"

template<typename REAL, typename SIZE>
class MeshRS;

template<typename REAL, typename SIZE>
std::ostream& operator<<(std::ostream& out, const MeshRS<REAL,SIZE>& mesh);

template<typename REAL, typename SIZE>
class FaceRS {
public:
  explicit FaceRS() { }
  explicit FaceRS(SIZE v0, SIZE v1, SIZE v2) { v[0] = v0; v[1] = v1; v[2] = v2; }
  SIZE  &operator[](SIZE i) { return v[i]; }
  SIZE   operator[](SIZE i) const { return v[i]; }

  Vector3R<REAL> nml; // face normal
private:
  SIZE v[3];
};

template<typename REAL>
class VertexR
{
public:
  explicit VertexR() : pos(0.0f,0.0f,0.0f), nml(0.0f,0.0f,0.0f) { }
  Vector3R<REAL> pos; // vertex position
  Vector3R<REAL> nml; // vertex normal
};

template<typename REAL>
using VertexVecR = std::vector< VertexR<REAL> >;

template<typename REAL, typename SIZE>
using FaceVecRS = std::vector< FaceRS<REAL, SIZE> >;

// defaults
typedef VertexVecR<double> VertexVec;
typedef FaceVecRS<double, unsigned int> FaceVec;

// A triangular mesh.
template<typename REAL, typename SIZE>
class MeshRS : public Primitive 
{
public:
  explicit MeshRS(const aiMesh *mesh);
  ~MeshRS();

  inline AlignedBox3f &compute_bbox();
  void compute_face_normals();

  const VertexVecR<REAL>      &get_verts() const { return m_verts; }
  const FaceVecRS<REAL, SIZE> &get_faces() const { return m_faces; }

  bool is_mesh() const { return true; }

  inline void transform_in_place(const AffineCompact3f &trans);

  friend std::ostream& operator<< <>(std::ostream& out, const MeshRS<REAL,SIZE>& mesh);

protected:
  VertexVecR<REAL>      m_verts;
  FaceVecRS<REAL, SIZE> m_faces;
};

template<typename REAL, typename SIZE>
using MeshPtrRS = boost::shared_ptr< MeshRS<REAL, SIZE> >;

// defaults
typedef MeshRS<double, unsigned int> Mesh;
typedef MeshPtrRS<double, unsigned int> MeshPtr;

#endif // MESH_H
