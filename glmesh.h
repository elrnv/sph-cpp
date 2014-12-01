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
      Mesh &mesh,
      bool dynamic,
      MaterialManager &matman,
      UniformBuffer &ubo,
      ShaderManager &shaderman);
  ~GLMesh();

  inline bool is_mesh()    const { return true; }

  inline Size get_num_indices()  const { return m_mesh.get_faces().size()*3; }
  inline Size get_num_vertices() const { return m_mesh.get_verts().size();   }

  // PRE: assume that point have been sorted
  virtual void update_glbuf_withsort(const AffineCompact3f &mvtrans, 
                                     const AffineCompact3f &nmltrans);
  virtual void update_glbuf_nosort();
  virtual void update_shader(ShaderManager::ShaderType type,
                             ShaderManager &shaderman);
  
  inline void print() const { std::cerr << m_mesh << std::endl; }

private: // helper functions
  void update_glbuf(const Matrix3XT<GLfloat> &vispos,
                    const Matrix3XT<GLfloat> &visnml);

private:
  Mesh &m_mesh; // reference to the native mesh object
};

typedef boost::shared_ptr< GLMesh > GLMeshPtr;

#endif // GLMESH_H
