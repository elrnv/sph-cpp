#include <QDebug>
#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLBuffer>
#include <limits>
#include <assimp/scene.h>
#include "mesh.h"

// Mesh stuff
template<typename REAL, typename SIZE>
MeshRS<REAL,SIZE>::MeshRS(const aiMesh *mesh)
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
  for ( auto &f : m_faces )
  {
    for (unsigned char i = 0; i < 3; ++i)
    {
      Vector3d e1 = (m_verts[f[(i+1)%3]].pos - m_verts[f[i]].pos).normalized();
      Vector3d e2 = (m_verts[f[(i+2)%3]].pos - m_verts[f[i]].pos).normalized();
      m_verts[f[i]].nml += e1.cross(e2); // compute the sum
    }
  }

  for ( auto &v : m_verts )
    v.nml.normalize();                  // take the average
}

template<typename REAL, typename SIZE>
MeshRS<REAL,SIZE>::~MeshRS()
{
}

template<typename REAL, typename SIZE>
void MeshRS<REAL,SIZE>::compute_face_normals()
{
  for ( auto &f : m_faces )
  {
    Vector3d e1 = m_verts[f[1]].pos - m_verts[f[0]].pos;
    Vector3d e2 = m_verts[f[2]].pos - m_verts[f[1]].pos;
    f.nml = e1.cross(e2).normalized();
  }
}

template<typename REAL, typename SIZE>
AlignedBox3f &MeshRS<REAL,SIZE>::compute_bbox()
{
  m_bbox.setEmpty();
  for ( auto &v : m_verts )
    m_bbox.extend(v.pos.template cast<float>());
  return m_bbox;
}

template<typename REAL, typename SIZE>
void MeshRS<REAL,SIZE>::transform_in_place(const AffineCompact3f &trans)
{
  // convert positions to canonical homogenized coordinates
  for ( auto &v : m_verts )
  {
    Vector3R<REAL> vec( trans.linear().template cast<REAL>() * v.pos );
    Translation3f T( trans.translation() );
    v.pos = vec + Vector3R<REAL>(T.x(), T.y(), T.z());
    v.nml = trans.inverse().linear().transpose().template cast<REAL>() * v.nml;
  }
  for ( auto &f : m_faces )
  {
    f.nml = trans.inverse().linear().transpose().template cast<REAL>() * f.nml;
  }
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

// defaults
template class MeshRS<double, unsigned int>;

template std::ostream& operator<<(std::ostream& out, const Mesh& mesh);
