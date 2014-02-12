#include <QDebug>
#include <limits>
#include <assimp/scene.h>
#include "mesh.h"

Mesh::Mesh(const aiMesh *mesh, bool need_bbox)
{
  // First copy all the vertices
  if (!mesh->HasPositions())
    return;

  aiVector3D *verts = mesh->mVertices;
  unsigned int num_verts = mesh->mNumVertices;
  m_verts.resize(num_verts);

  VertexVec::iterator v_it = m_verts.begin(); // member face iterator
  for (unsigned int i = 0; i < num_verts; ++i, ++v_it)
    v_it->pos << verts[i].x, verts[i].y, verts[i].z;

  if (need_bbox)
    compute_bbox(); // compute the bounded box if requested

  if (mesh->HasNormals())
  {
    aiVector3D *normals = mesh->mNormals;
    v_it = m_verts.begin();
    for (unsigned int i = 0; i < num_verts; ++i, ++v_it)
      v_it->nml << normals[i].x, normals[i].y, normals[i].z;
  }

  if (!mesh->HasFaces())
    return;

  aiFace *faces = mesh->mFaces;
  unsigned int num_faces = mesh->mNumFaces;
  m_faces.resize(num_faces);
  
  // Assume all faces are already triangles as returned by assimp
  FaceVec::iterator f_it = m_faces.begin(); // member face iterator
  for (unsigned int i = 0; i < num_faces; ++i, ++f_it)
  {
    if (faces[i].mNumIndices < 3)
      continue; // ignore line and point primitives
    (*f_it)[0] = SIZE_IDX(faces[i].mIndices[0]);
    (*f_it)[1] = SIZE_IDX(faces[i].mIndices[1]);
    (*f_it)[2] = SIZE_IDX(faces[i].mIndices[2]);
  }

  if (mesh->HasNormals())
    return; // Already collected the normals

  // Each vertex normal is the average of the neighbouring face normals
  for (f_it = m_faces.begin(); f_it != m_faces.end(); ++f_it)
  {
    Face &face = *f_it;
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

void Mesh::compute_face_normals()
{
  for (FaceVec::iterator f_it = m_faces.begin(); f_it != m_faces.end(); ++f_it)
  {
    Face &face = *f_it;
    Vector3d e1 = m_verts[face[1]].pos - m_verts[face[0]].pos;
    Vector3d e2 = m_verts[face[2]].pos - m_verts[face[1]].pos;
    face.nml = e1.cross(e2).normalized();
  }
}

void Mesh::compute_bbox()
{
  m_bbox.setEmpty();
  VertexVec::iterator v_it;
  for (v_it= m_verts.begin(); v_it != m_verts.end(); ++v_it)
    m_bbox.extend(v_it->pos.cast<float>());
}

std::ostream& operator<<(std::ostream& out, const Mesh& mesh)
{
  out << "mesh({ ";
  VertexVec::const_iterator v_it;
  for (v_it = mesh.m_verts.begin(); v_it != mesh.m_verts.end(); ++v_it)
  {
    if (v_it != mesh.m_verts.begin()) out << ",\n       ";
    out << v_it->pos[0] << " " << v_it->pos[1] << " " << v_it->pos[2];
  }
  out << "},\n\n     { ";

  FaceVec::const_iterator f_it;
  for (f_it = mesh.m_faces.begin(); f_it != mesh.m_faces.end(); ++f_it)
  {
    if (f_it != mesh.m_faces.begin()) out << ",\n       ";
    out << "[ " << (*f_it)[0] << " " << (*f_it)[1] << " " << (*f_it)[2] << " ]";
  }
  out << "});\n";
  return out;
}

