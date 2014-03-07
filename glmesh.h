#ifndef MESH_H
#define MESH_H

#include <vector>
#include <iostream>
#include <assimp/mesh.h>
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

template<typename SIZE>
using FaceVecS = FaceVecRS<double, SIZE>;

// defaults
typedef VertexVecR<double> VertexVec;
typedef FaceVecRS<double, unsigned int> FaceVec;

// A triangular mesh.
template<typename REAL, typename SIZE>
class MeshRS : public Primitive 
{
public:
  explicit MeshRS(const aiMesh *mesh, bool compute_bbox = true);
  ~MeshRS();

  void compute_bbox();
  void compute_face_normals();

  const VertexVecR<REAL>      &get_verts() const { return m_verts; }
  const FaceVecRS<REAL, SIZE> &get_faces() const { return m_faces; }

  bool is_mesh() const { return true; }

  friend std::ostream& operator<< <>(std::ostream& out, const MeshRS<REAL,SIZE>& mesh);

protected:
  VertexVecR<REAL>      m_verts;
  FaceVecRS<REAL, SIZE> m_faces;
};

template<typename SIZE>
using MeshS = MeshRS<double, SIZE>;

// A triangular mesh representation for OpenGL applications
template<typename SIZE>
class GLMeshS : public GLPrimitiveS<SIZE>
{
public:
  explicit GLMeshS(
      MeshS<SIZE> *mesh,
      const Material *mat,
      UniformBuffer &ubo,
      ShaderManager &shaderman);
  ~GLMeshS();

  MeshS<SIZE> *get_mesh() { return m_mesh; }

  bool is_mesh() const { return true; }

  SIZE get_num_indices()  const { return m_mesh->get_faces().size()*3; }
  SIZE get_num_vertices() const { return m_mesh->get_verts().size();   }

  void update_data();
  void update_glbuf();
  void update_shader(ShaderManager::ShaderType type);
  
  void print() const { std::cerr << *m_mesh << std::endl; }

protected:
  MeshS<SIZE> *m_mesh; // reference to the native mesh object

  // intermediate buffers between mesh and glbuffers
  std::vector<GLfloat> m_vertices; 
  std::vector<GLfloat> m_normals;
  bool                 m_insync;
};

// defaults
typedef MeshRS<double, unsigned int> Mesh;
typedef GLMeshS<unsigned int> GLMesh;

#endif // MESH_H
