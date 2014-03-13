#include <QDebug>
#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLBuffer>
#include <limits>
#include <assimp/scene.h>
#include "glmesh.h"

// GLMesh stuff
template<typename REAL, typename SIZE>
GLMeshRS<REAL,SIZE>::GLMeshRS(
    MeshRS<REAL,SIZE> *mesh,
    const Material *mat,
    UniformBuffer &ubo,
    ShaderManager &shaderman)
  : GLPrimitiveS<SIZE>(mat, ubo, shaderman)
  , m_mesh(mesh)
  , m_insync(true)
{
  // collect vertex attributes
  const VertexVec &verts = mesh->get_verts();

  GLfloat vertices[verts.size()*3];
  GLfloat normals[verts.size()*3];

  VertexVec::const_iterator v_it = verts.cbegin();
  int i = 0;
  for ( ; v_it != verts.cend(); ++v_it)
    for (short j = 0; j < 3; ++j)
    {
      vertices[i] = GLfloat(v_it->pos[j]);
      normals[i++]  = GLfloat(v_it->nml[j]);
    }

  // collect indices
  const FaceVecRS<REAL, SIZE> &faces = mesh->get_faces();

  GLuint indices[faces.size()*3];

  typename FaceVecRS<REAL, SIZE>::const_iterator f_it = faces.begin();
  i = 0;
  for ( ; f_it != faces.end(); ++f_it)
    for (short j = 0; j < 3; ++j)
    {
      indices[i] = GLuint((*f_it)[j]);
      ++i;
    }

  this->m_vao.create();
  this->m_vao.bind();

  this->m_pos.create();
  this->m_pos.setUsagePattern( QOpenGLBuffer::StaticDraw );
  this->m_pos.bind();
  this->m_pos.allocate( vertices, sizeof( vertices ) );

  this->m_nml.create();
  this->m_nml.setUsagePattern( QOpenGLBuffer::StaticDraw );
  this->m_nml.bind();
  this->m_nml.allocate( normals, sizeof( normals ) );

  this->m_idx.create();
  this->m_idx.setUsagePattern( QOpenGLBuffer::StaticDraw );
  this->m_idx.bind();
  this->m_idx.allocate( indices, sizeof( indices ) );

  this->m_vao.release();
  update_shader(ShaderManager::PHONG);
}

template<typename REAL, typename SIZE>
GLMeshRS<REAL,SIZE>::~GLMeshRS()
{
}

template<typename REAL, typename SIZE>
void GLMeshRS<REAL, SIZE>::update_data()
{
  std::lock_guard<std::mutex> guard(this->m_lock);

  if (!m_insync)
    return;

  // collect vertex attributes
  const VertexVec &verts = m_mesh->get_verts();

  m_vertices.resize(verts.size()*3);
  m_normals.resize(verts.size()*3);

  int i = 0;
  for ( auto &v : verts )
  {
    for (short j = 0; j < 3; ++j)
    {
      m_vertices[i] = GLfloat(v.pos[j]);
      m_normals[i++] = GLfloat(v.nml[j]);
    }
  }

  m_insync = false;
}

template<typename REAL, typename SIZE>
void GLMeshRS<REAL, SIZE>::update_glbuf()
{
  std::lock_guard<std::mutex> guard(this->m_lock);
  
  if (m_insync) //TODO: test with atomic m_insync
    return;

  this->m_vao.bind();
  this->m_pos.bind();
  this->m_pos.write( 0, m_vertices.data(), sizeof( m_vertices ) );

  this->m_nml.bind();
  this->m_nml.write( 0, m_normals.data(), sizeof( m_normals ) );
  this->m_vao.release();

  m_insync = true;
}

template<typename REAL, typename SIZE>
void GLMeshRS<REAL, SIZE>::update_shader(ShaderManager::ShaderType type)
{
  if (this->m_prog)
  {
    this->m_prog->disableAttributeArray( "pos" );
    this->m_prog->disableAttributeArray( "nml" );
  }

  if (type == ShaderManager::PHONG)
    this->m_prog = this->m_shaderman.get_phong_shader();
  else if (type == ShaderManager::PARTICLE)
    this->m_prog = this->m_shaderman.get_particle_shader();
  else
    this->m_prog = this->m_shaderman.get_wireframe_shader();

  this->m_vao.bind();

  this->m_pos.bind();
  this->m_prog->enableAttributeArray( "pos" );
  this->m_prog->setAttributeBuffer( "pos", GL_FLOAT, 0, 3 );

  this->m_nml.bind();
  this->m_prog->enableAttributeArray( "nml" );
  this->m_prog->setAttributeBuffer( "nml", GL_FLOAT, 0, 3 );

  this->m_vao.release();

  this->m_ubo.bindToProg(this->m_prog->programId(), "Globals");
}

// defaults
template class GLMeshRS<double, unsigned int>;
