#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <assimp/mesh.h>
#include <assimp/scene.h>
#include "eigen.h"
#include "primitive.h"
#include "material.h"
#include "dynparams.h"

class Scene
{
private:
  class Node 
  {
  public:
    Node(const std::string& name);
    Node(const aiNode *node); // construct a scene node from an assimp node
    ~Node();

    const std::string &get_name() const { return m_name; } // for debugging
    const AffineCompact3f& get_trans() const { return m_trans; }

    // rotation in degrees
    void rotate(float angle, const Vector3f &axis);

    void scale(float amount); // uniform scale
    void scale(const Vector3f& amount); // non-uniform scale
    void translate(const Vector3f& amount);

    // compute the parent transformations into the children, reset m_trans
    // this operation invalidates the bounding box
    void flatten();

    void print(int depth = 0) const;

    bool is_geometry() const { return false; }

    unsigned int num_primitives() const;

    AlignedBox3f &compute_bbox();
    AlignedBox3f &get_bbox() { return m_bbox; }
    void cube_bbox();

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
    std::string     m_name;      // human readable name
    AffineCompact3f m_trans;     // transformations
    AlignedBox3f    m_bbox;      // bounding box
    Index           m_prim_idx;  // optional index to the represented primitive
    NodeType        m_node_type; // what kind of node is this

    // child data (indices into the main node list)
    Index m_first_child;
    Size  m_num_children;
  };

  typedef std::vector<Node> NodeVec;

public:
  Scene(const aiScene *scene); // construct a scene from an assimp scene
  ~Scene();

  void add_node(const aiNode *ainode);

private:
  NodeVec m_nodes;
};


class GeometryNode : public Node 
{
public:
  // Create a geo node with existing primitive
  // Ownership of the primitive is transferred to this node
  GeometryNode(const std::string& name, PrimitivePtr primitive);

  // Construct a geometry node and mesh from given assimp mesh
  GeometryNode( const aiMesh *mesh, const aiMaterial *mat );

  // copy constructor
  GeometryNode(const GeometryNode &orig);

  // virtual allocation routine to avoid using dynamic_casts
  virtual GeometryNode *clone() const { return new GeometryNode(*this); }

  virtual ~GeometryNode();

  virtual bool is_geometry() const { return true; }
  virtual bool is_dynamic()  const { return false; }

  virtual PrimitivePtr     get_primitive() const { return m_primitive; }
          MaterialConstPtr get_material()  const { return m_material; }

  // compute the transformations into the primitives
  // this operation invalidates the bounding box
  virtual void flatten();

  void print(int depth = 0) const;

  AlignedBox3f &compute_bbox();
  void cube_bbox();

protected:
  PrimitivePtr      m_primitive;
  MaterialConstPtr  m_material;
};


class FluidNode : public GeometryNode 
{
public:
  // Construct a geometry node and mesh from given assimp mesh
  FluidNode(
      const aiMesh *mesh,
      const aiMaterial *mat,
      DynParamsPtr dyn_params);

  FluidNode(const FluidNode &orig); // copy constructor
  virtual FluidNode *clone() const { return new FluidNode(*this); }
  virtual ~FluidNode();

  virtual bool is_dynamic() const { return true; }

  // using covariant return type since primitive is encapsulated by m_fluid
  virtual FluidPtr get_primitive() const { return m_fluid; }

protected:
  FluidPtr m_fluid;
}; // class FluidNode

#if 0
class RigidNode : public GeometryNode 
{
public:
  // Construct a geometry node and mesh from given assimp mesh
  RigidNode(
      const aiMesh *mesh,
      const aiMaterial *mat,
      DynParamsPtr dyn_params);

  RigidNode(const RigidNode &orig); // copy constructor
  virtual RigidNode *clone() const { return new RigidNode(*this); }
  virtual ~RigidNode();

  virtual bool is_dynamic() const { return true; }

  virtual RigidPtr get_primitive() const { return m_rigid; }

protected:
  RigidPtr m_rigid;
}; // class RigidNode
#endif

#endif // SCENE_H
