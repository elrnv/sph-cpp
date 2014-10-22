#ifndef GLMESH_H
#define GLMESH_H

#include <vector>
#include <iostream>
#include <assimp/mesh.h>
#include <atomic>
#include "mesh.h"
#include "glprimitive.h"

// A triangular mesh representation for OpenGL applications

class GLMesh : public GLPrimitive
{
public:
  explicit GLMesh(
      MeshPtr mesh,
      MaterialConstPtr mat,
      UniformBuffer &ubo,
      ShaderManager &shaderman);
  ~GLMesh();

  inline MeshPtr get_mesh() { return m_mesh; }

  inline bool is_mesh()    const { return true; }

  inline Size get_num_indices()  const { return m_mesh->get_faces().size()*3; }
  inline Size get_num_vertices() const { return m_mesh->get_verts().size();   }

  void update_data();

  inline Vector3f get_closest_pt() const { return Vector3f(m_vertices.col(0)); }
  void sort_by_depth(const AffineCompact3f &mvtrans);
  void update_glbuf();
  void update_shader(ShaderManager::ShaderType type);
  
  inline void print() const { std::cerr << *m_mesh << std::endl; }

private:
  MeshPtr m_mesh; // reference to the native mesh object

  // intermediate buffers between mesh and glbuffers
  Matrix3XR<GLfloat> m_vertices; 
  Matrix3XR<GLfloat> m_normals;
  std::atomic<bool>  m_insync;
};

#endif // GLMESH_H
