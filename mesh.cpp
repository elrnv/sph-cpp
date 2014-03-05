#include <QDebug>
#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLBuffer>
#include <limits>
#include <assimp/scene.h>
#include "mesh.h"

// Mesh stuff

template<typename REAL, typename SIZE>
MeshRS<REAL,SIZE>::MeshRS(const aiMesh *mesh, bool need_bbox)
{
  // First copy all the vertices
  if (!mesh->HasPositions())
    return;

  aiVector3D *verts = mesh->mVertices;
  SIZE num_verts = mesh->mNumVertices;
  m_verts.resize(num_verts);

  typename VertexVecR<REAL>::iterator v_it = m_verts.begin(); // member face iterator
  for (SIZE i = 0; i < num_verts; ++i, ++v_it)
    v_it->pos << verts[i].x, verts[i].y, verts[i].z;

  if (need_bbox)
    compute_bbox(); // compute the bounded box if requested

  if (mesh->HasNormals())
  {
    aiVector3D *normals = mesh->mNormals;
    v_it = m_verts.begin();
    for (SIZE i = 0; i < num_verts; ++i, ++v_it)
      v_it->nml << normals[i].x, normals[i].y, normals[i].z;
  }

  if (!mesh->HasFaces())
    return;

  aiFace *faces = mesh->mFaces;
  SIZE num_faces = mesh->mNumFaces;
  m_faces.reserve(num_faces);
  
  // Assume all faces are already triangles as returned by assimp
  for (SIZE i = 0; i < num_faces; ++i)
  {
    if (faces[i].mNumIndices < 3)
      continue; // ignore line and point primitives
    m_faces.push_back(
        FaceRS<REAL, SIZE>(
          faces[i].mIndices[0],
          faces[i].mIndices[1],
          faces[i].mIndices[2]));
  }

  if (mesh->HasNormals())
    return; // Already collected the normals

  // Each vertex normal is the average of the neighbouring face normals
  typename FaceVecRS<REAL,SIZE>::iterator f_it = m_faces.begin(); // member face iterator
  for (f_it = m_faces.begin(); f_it != m_faces.end(); ++f_it)
  {
    FaceRS<REAL,SIZE> &face = *f_it;
    for (unsigned char i = 0; i < 3; ++i)
    {
      Vector3d e1 = (m_verts[face[(i+1)%3]].pos - m_verts[face[i]].pos).normalized();
      Vector3d e2 = (m_verts[face[(i+2)%3]].pos - m_verts[face[i]].pos).normalized();
      m_verts[(*f_it)[i]].nml += e1.cross(e2); // compute the sum
    }
  }

  for (v_it = m_verts.begin(); v_it != m_verts.end(); ++v_it)
    v_it->nml.normalize();                  // take the average
}

template<typename REAL, typename SIZE>
MeshRS<REAL,SIZE>::~MeshRS()
{
}

template<typename REAL, typename SIZE>
void MeshRS<REAL,SIZE>::compute_face_normals()
{
  typename FaceVecRS<REAL,SIZE>::iterator f_it = m_faces.begin();
  for ( ; f_it != m_faces.end(); ++f_it)
  {
    FaceRS<REAL,SIZE> &face = *f_it;
    Vector3d e1 = m_verts[face[1]].pos - m_verts[face[0]].pos;
    Vector3d e2 = m_verts[face[2]].pos - m_verts[face[1]].pos;
    face.nml = e1.cross(e2).normalized();
  }
}

template<typename REAL, typename SIZE>
void MeshRS<REAL,SIZE>::compute_bbox()
{
  m_bbox.setEmpty();
  typename VertexVecR<REAL>::iterator v_it;
  for (v_it = m_verts.begin(); v_it != m_verts.end(); ++v_it)
    m_bbox.extend(v_it->pos.template cast<float>());
}

template<typename REAL, typename SIZE>
std::ostream& operator<<(std::ostream& out, const MeshRS<REAL,SIZE>& mesh)
{
  out << "mesh({ ";
  typename VertexVecR<REAL>::const_iterator v_it;
  for (v_it = mesh.m_verts.begin(); v_it != mesh.m_verts.end(); ++v_it)
  {
    if (v_it != mesh.m_verts.begin()) out << ",\n       ";
    out << v_it->pos[0] << " " << v_it->pos[1] << " " << v_it->pos[2];
  }
  out << "},\n\n     { ";

  typename FaceVecRS<REAL,SIZE>::const_iterator f_it;
  for (f_it = mesh.m_faces.begin(); f_it != mesh.m_faces.end(); ++f_it)
  {
    if (f_it != mesh.m_faces.begin()) out << ",\n       ";
    out << "[ " << (*f_it)[0] << " " << (*f_it)[1] << " " << (*f_it)[2] << " ]";
  }
  out << "});\n";
  return out;
}


// GLMesh stuff

template<typename SIZE>
GLMeshS<SIZE>::GLMeshS(
    MeshS<SIZE> *mesh,
    const PhongMaterial *mat,
    UniformBuffer &ubo,
    ShaderManager &shaderman)
  : GLPrimitiveS<SIZE>(mat, ubo, shaderman)
  , m_mesh(mesh)
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
  const FaceVecS<SIZE> &faces = mesh->get_faces();

  GLuint indices[faces.size()*3];

  typename FaceVecS<SIZE>::const_iterator f_it = faces.begin();
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

  update_shader(ShaderManager::PHONG);
}

template<typename SIZE>
GLMeshS<SIZE>::~GLMeshS()
{
}


template<typename SIZE>
void GLMeshS<SIZE>::update_shader(ShaderManager::ShaderType type)
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

  this->m_ubo.bindToProg(this->m_prog->programId(), "Globals");
}


template class MeshRS<double, unsigned int>;
template class GLMeshS<unsigned int>;

template std::ostream& operator<<(std::ostream& out, const Mesh& mesh);
