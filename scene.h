#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <assimp/mesh.h>
#include <assimp/scene.h>
#include "eigen.h"
#include "primitive.h"
#include "material.h"

class SceneNode;
typedef std::vector<SceneNode*> NodeList;

class SceneNode 
{
public:

  SceneNode(const std::string& name);

  // construct a scene node from an assimp node
  SceneNode(const aiNode *node);

  // copy constructor for deep copy
  SceneNode(const SceneNode &orig);

  // virtual allocation routine to avoid using dynamic_casts
  virtual SceneNode *clone() { return new SceneNode(*this); }

  virtual ~SceneNode();

  const std::string &get_name() const { return m_name; } // for debugging

  const AffineCompact3d& get_trans() const { return m_trans; }

  // ownership of child is implicitly passed to this node
  void add_child(SceneNode* child)
  {
    m_children.push_back(child);
  }

  void rotate(char axis, double angle);
  void scale(const Vector3d& amount);
  void translate(const Vector3d& amount);

  // precompute the parent transformations into the children 
  void flatten();

  virtual void print(int depth = 0) const;

  NodeList::const_iterator begin() { return m_children.begin(); }
  NodeList::const_iterator end()   { return m_children.end(); }

  virtual bool is_geometry() const { return false; }

protected:
  // Useful for picking
  std::string m_name;

  // Transformations
  AffineCompact3d m_trans;

  // Hierarchy
  NodeList m_children;
};


class GeometryNode : public SceneNode 
{
public:
  // Create a geo node with existing primitive
  // Ownership of the primitive is transferred to this node
  GeometryNode(const std::string& name, Primitive* primitive);

  // Construct a geometry node and mesh from given assimp mesh
  GeometryNode(const aiMesh *mesh);

  // copy constructor
  GeometryNode(const GeometryNode &orig);

  // virtual allocation routine to avoid using dynamic_casts
  virtual GeometryNode *clone() { return new GeometryNode(*this); }

  virtual ~GeometryNode();

  virtual bool is_geometry() const { return true; }

  const Primitive* get_primitive() const { return m_primitive; }
  const Material*  get_material()  const { return m_material; }

  // The following two routines return overwritten memebers
  Material *set_material(Material* material)
  {
    Material *temp = m_material;
    m_material = material;
    return temp;
  }

  Primitive *set_primitive(Primitive* primitive)
  {
    Primitive *temp = m_primitive;
    m_primitive = primitive;
    return temp;
  }

  void print(int depth = 0) const;

protected:
  Primitive *m_primitive;
  Material  *m_material;

};

#endif // SCENE_H
