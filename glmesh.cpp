#include <QDebug>
#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLBuffer>
#include <limits>
#include <assimp/scene.h>
#include "glmesh.h"

// GLMesh stuff
template<typename REAL, typename SIZE>
GLMeshRS<REAL,SIZE>::GLMeshRS(
    MeshPtrRS<REAL,SIZE> mesh,
    MaterialConstPtr mat,
    UniformBuffer &ubo,
    ShaderManager &shaderman)
  : GLPrimitiveS<SIZE>(mat, ubo, shaderman)
  , m_mesh(mesh)
  , m_vertices(3, get_num_vertices())
  , m_normals(3, get_num_vertices())
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

  //update_data(); // TODO: uncomment when sort_by_depth is done

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

  int i = 0;
  for ( auto &v : verts )
  {
    m_vertices.col(i) = v.pos.template cast<GLfloat>();
    m_normals.col(i)  = v.nml.template cast<GLfloat>();
    i += 1;
  }

  m_insync = false;
}

template<typename REAL, typename SIZE>
void GLMeshRS<REAL, SIZE>::sort_by_depth(const AffineCompact3f &mvtrans)
{
  return;
  // TODO: implement this
  std::lock_guard<std::mutex> guard(this->m_lock);

  // Sort all vertices by the z value
  Matrix3XR<GLfloat> mvpos = (mvtrans * m_vertices).eval(); // TODO: mem alloc expensive?
  SIZE num_verts = get_num_vertices();
  VectorXT<SIZE> perm_vec(num_verts);
  for (SIZE i = 0; i < num_verts; ++i)
    perm_vec[i] = i;

  SIZE *perm_data = perm_vec.data();

  std::sort(perm_data, perm_data + num_verts,
      [mvpos](SIZE i, SIZE j) { return mvpos.col(i)[2] < mvpos.col(j)[2]; });

  PermutationMatrix<Dynamic, Dynamic, SIZE> perm_mat(perm_vec);
  m_vertices = m_vertices * perm_mat;
  m_normals  = m_normals * perm_mat;

  for (SIZE i = 0; i < num_verts; ++i)
    std::cerr << m_vertices.col(i)[0] << " " << m_vertices.col(i)[1] << " " << m_vertices.col(i)[2] << std::endl;

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
    this->m_prog->bind();
    this->m_prog->disableAttributeArray( "pos" );
    this->m_prog->disableAttributeArray( "nml" );
  }

  if (type == ShaderManager::PHONG)
    this->m_prog = this->m_shaderman.get_phong_shader();
  else if (type == ShaderManager::PARTICLE)
    this->m_prog = this->m_shaderman.get_particle_shader();
  else
    this->m_prog = this->m_shaderman.get_wireframe_shader();

  this->m_prog->bind();
  this->m_vao.bind();

  this->m_pos.bind();
  this->m_prog->enableAttributeArray( "pos" );
  this->m_prog->setAttributeBuffer( "pos", GL_FLOAT, 0, 3 );

  this->m_nml.bind();
  this->m_prog->enableAttributeArray( "nml" );
  this->m_prog->setAttributeBuffer( "nml", GL_FLOAT, 0, 3 );

  this->m_vao.release();

  this->m_ubo.bindToProg(this->m_prog->programId(), "Globals");
  this->m_prog->release();
}

// defaults
template class GLMeshRS<double, unsigned int>;
