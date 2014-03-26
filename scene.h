#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <assimp/mesh.h>
#include <assimp/scene.h>
#include "eigen.h"
#include "primitive.h"
#include "material.h"
#include "util.h"

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
  virtual SceneNode *clone() const { return new SceneNode(*this); }

  virtual ~SceneNode();

  const std::string &get_name() const { return m_name; } // for debugging

  const AffineCompact3f& get_trans() const { return m_trans; }

  // ownership of child is implicitly passed to this node
  void add_child(SceneNode* child)
  {
    m_children.push_back(child);
  }

  // rotation in degrees
  void rotate(float angle, const Vector3f &axis);

  void scale(float amount); // uniform scale
  void scale(const Vector3f& amount); // non-uniform scale
  void translate(const Vector3f& amount);

  // compute the parent transformations into the children, reset m_trans
  // this operation invalidates the bounding box
  virtual void flatten();

  virtual void print(int depth = 0) const;

  const NodeList &get_children() const { return m_children; }

  virtual bool is_geometry() const { return false; }

  unsigned int num_primitives() const;

  virtual AlignedBox3f &compute_bbox();
  AlignedBox3f &get_bbox() { return m_bbox; }
  virtual void cube_bbox();

  // translate and scale model to fit in a 2x2x2 box centered at the origin
  void normalize_model();

  // same as above but extend the bounding box according to given vectors
  // in the form Vector2f( positive extension, negative extension )
  void normalize_model(
      const Vector2f &ext_x,
      const Vector2f &ext_y,
      const Vector2f &ext_z);

  // same as normalize_model(void) but extending the bottom left floor corner
  // and top right ceil corner by the given vectors
  void normalize_model(
      const Vector3f &ext_blf,
      const Vector3f &ext_trc);

protected:
  // Useful for picking
  std::string m_name;

  // Transformations
  AffineCompact3f m_trans;

  // Hierarchy
  NodeList m_children;

  AlignedBox3f m_bbox;
};


class GeometryNode : public SceneNode 
{
public:
  // Create a geo node with existing primitive
  // Ownership of the primitive is transferred to this node
  GeometryNode(const std::string& name, Primitive* primitive);

  // Construct a geometry node and mesh from given assimp mesh
  GeometryNode(
      const aiMesh *mesh,
      const aiMaterial *mat,
      const DynParams &dyn_params);

  // copy constructor
  GeometryNode(const GeometryNode &orig);

  // virtual allocation routine to avoid using dynamic_casts
  virtual GeometryNode *clone() const { return new GeometryNode(*this); }

  virtual ~GeometryNode();

  virtual bool is_geometry() const { return true; }

        Primitive* get_primitive() const { return m_primitive; }
  const Material*  get_material()  const { return m_material; }

  // The following two routines return overwritten memebers
  const Material *set_material(const Material* material)
  {
    const Material *temp = m_material;
    m_material = material;
    return temp;
  }

  Primitive *set_primitive(Primitive* primitive)
  {
    Primitive *temp = m_primitive;
    m_primitive = primitive;
    return temp;
  }

  // compute the transformations into the primitives
  // this operation invalidates the bounding box
  virtual void flatten();

  void print(int depth = 0) const;

  AlignedBox3f &compute_bbox();
  void cube_bbox();

protected:
  Primitive *m_primitive;
  const Material *m_material;

};

#endif // SCENE_H
