#ifndef MESH_H
#define MESH_H

#include <vector>
#include <assimp/mesh.h>
#include "primitive.h"
#include "eigen.h"

typedef unsigned short SIZE_IDX;

class Face {
public:
  explicit Face() { }
  explicit Face(SIZE_IDX v0, SIZE_IDX v1, SIZE_IDX v2) { v[0] = v0; v[1] = v1; v[2] = v2; }
  SIZE_IDX  &operator[](SIZE_IDX i) { return v[i]; }
  SIZE_IDX   operator[](SIZE_IDX i) const { return v[i]; }

  Vector3d nml; // face normal
private:
  SIZE_IDX v[3];
};

struct Vertex
{
  explicit Vertex() : pos(0.0f,0.0f,0.0f), nml(0.0f,0.0f,0.0f) { }
  Vector3d pos; // vertex position
  Vector3d nml; // vertex normal
};

typedef std::vector<Vertex> VertexVec;
typedef std::vector<Face>   FaceVec;

// A triangular mesh.
class Mesh : public Primitive 
{
public:
  explicit Mesh(const aiMesh *mesh, bool compute_bbox = true);

  void compute_bbox();
  void compute_face_normals();

  const VertexVec &get_verts() const { return m_verts; }
  const FaceVec   &get_faces() const { return m_faces; }

  bool is_mesh() const { return true; }

protected:
  VertexVec    m_verts;
  FaceVec      m_faces;
  AlignedBox3f m_bbox;

  friend std::ostream& operator<<(std::ostream& out, const Mesh& mesh);
};

#endif // MESH_H
