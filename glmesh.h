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

  inline Vector3f get_closest_pt() const { return Vector3f(m_vertices.col(0)); }
  void update_glbuf_withsort(const AffineCompact3f &mvtrans);
  void update_glbuf_nosort();
  void update_shader(ShaderManager::ShaderType type);
  
  inline void print() const { std::cerr << *m_mesh << std::endl; }

private:
  MeshPtr m_mesh; // reference to the native mesh object
};

#endif // GLMESH_H
