#ifndef GLMESH_H
#define GLMESH_H

#include <vector>
#include <iostream>
#include <assimp/mesh.h>
#include "mesh.h"
#include "glprimitive.h"

// A triangular mesh representation for OpenGL applications
template<typename REAL, typename SIZE>
class GLMeshRS : public GLPrimitiveS<SIZE>
{
public:
  explicit GLMeshRS(
      MeshRS<REAL,SIZE> *mesh,
      const Material *mat,
      UniformBuffer &ubo,
      ShaderManager &shaderman);
  ~GLMeshRS();

  MeshRS<REAL,SIZE> *get_mesh() { return m_mesh; }

  bool is_mesh() const { return true; }

  SIZE get_num_indices()  const { return m_mesh->get_faces().size()*3; }
  SIZE get_num_vertices() const { return m_mesh->get_verts().size();   }

  void update_data();
  void update_glbuf();
  void update_shader(ShaderManager::ShaderType type);
  
  void print() const { std::cerr << *m_mesh << std::endl; }

protected:
  MeshRS<REAL,SIZE> *m_mesh; // reference to the native mesh object

  // intermediate buffers between mesh and glbuffers
  std::vector<GLfloat> m_vertices; 
  std::vector<GLfloat> m_normals;
  bool                 m_insync;
};

// defaults
typedef GLMeshRS<double, unsigned int> GLMesh;

#endif // GLMESH_H