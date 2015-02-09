#include <QDebug>
#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLBuffer>
#include <limits>
#include <assimp/scene.h>
#include "materialmanager.h"
#include "glmesh.h"

// GLMesh stuff

GLMesh::GLMesh(
    Mesh &mesh,
    bool dynamic,
    MaterialManager &matman,
    UniformBuffer &ubo,
    ShaderManager &shaderman)
  : GLPrimitive(matman[mesh.get_material_idx()], ubo)
  , m_mesh(mesh)
{
  (void) dynamic; // TODO: implement dynamic meshes
  // collect vertex attributes
  const VertexVec &verts = mesh.get_verts();

  GLfloat vertices[verts.size()*3];
  GLfloat normals[verts.size()*3];

  VertexVec::const_iterator v_it = verts.cbegin();
  int i = 0;
  for ( ; v_it != verts.cend(); ++v_it)
  {
    for (short j = 0; j < 3; ++j)
    {
      vertices[i] = GLfloat(v_it->pos[j]);
      normals[i++]  = GLfloat(v_it->nml[j]);
    }
  }

  // collect indices
  const FaceVec &faces = mesh.get_faces();

  GLuint indices[faces.size()*3];

  typename FaceVec::const_iterator f_it = faces.begin();
  i = 0;
  for ( ; f_it != faces.end(); ++f_it)
  {
    for (short j = 0; j < 3; ++j)
    {
      indices[i] = GLuint((*f_it)[j]);
      ++i;
    }
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

  update_shader(ShaderManager::PHONG, shaderman);
}


GLMesh::~GLMesh()
{
}

void
GLMesh::update_glbuf(const Matrix3XT<GLfloat> &vispos,
                     const Matrix3XT<GLfloat> &visnml)
{
  this->m_vao.bind();
  this->m_pos.bind();
  this->m_pos.write( 0, vispos.data(), sizeof( vispos ) );

  this->m_nml.bind();
  this->m_nml.write( 0, visnml.data(), sizeof( visnml ) );
  this->m_vao.release();
}

void 
GLMesh::update_glbuf_withsort(const AffineCompact3f &mvtrans,
                              const AffineCompact3f &nmlmvtrans)
{ // no actual sort needed here
  if (!m_mesh.is_staleposnml())
  {
    update_glbuf(mvtrans * m_mesh.get_vispos(),
                 nmlmvtrans * m_mesh.get_visnml());
    m_mesh.set_staleposnml(true);
  }
}

void
GLMesh::update_glbuf_nosort()
{
  if (!m_mesh.is_staleposnml())
  {
    update_glbuf(m_mesh.get_vispos(), m_mesh.get_visnml());
    m_mesh.set_staleposnml(true);
  }
}


void
GLMesh::update_shader(ShaderManager::ShaderType type,
                      ShaderManager &shaderman)
{
  if (this->m_prog)
  {
    this->m_prog->bind();
    this->m_prog->disableAttributeArray( "pos" );
    this->m_prog->disableAttributeArray( "nml" );
  }

  if (type == ShaderManager::PHONG)
    this->m_prog = shaderman.get_phong_shader();
  else if (type == ShaderManager::PARTICLE)
    this->m_prog = shaderman.get_particle_shader();
  else
    this->m_prog = shaderman.get_normals_shader();

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
