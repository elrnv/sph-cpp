#include <iostream>
#include <limits>
#include <Eigen/Dense>
#include "scene.h"
#include "mesh.h"

SceneNode::SceneNode(const std::string& name)
  : m_name(name)
{
  m_trans.setIdentity();
}

SceneNode::SceneNode(const aiNode *node)
  : m_name(node->mName.C_Str())
{
  const aiMatrix4x4 &trans = node->mTransformation;
  // verify that trans is affine
  if (trans.d1 != 0 || trans.d2 != 0 || trans.d3 != 0 || trans.d4 != 1)
    fprintf(stderr, "WARNING: Loaded file contains non-affine transforms!");

  Matrix4d mat;
  mat << trans.a1, trans.a2, trans.a3, trans.a4,
         trans.b1, trans.b2, trans.b3, trans.b4,
         trans.c1, trans.c2, trans.c3, trans.c4,
         trans.d1, trans.d2, trans.d3, trans.d4;
  m_trans = mat;
}

SceneNode::SceneNode(const SceneNode &orig)
  : m_name(orig.m_name)
  , m_trans(orig.m_trans)
{
  m_children.resize(orig.m_children.size());

  // deep copy of children
  for( NodeList::const_iterator it = orig.m_children.begin();
       it != orig.m_children.end();
       it++ )
  {
    SceneNode *new_node = (*it)->clone(); // dynamic allocation with new
    m_children.push_back(new_node);
  }
}

SceneNode::~SceneNode()
{
}

void SceneNode::rotate(char axis, double angle)
{
}

void SceneNode::scale(const Vector3d& amount)
{
}

void SceneNode::translate(const Vector3d& amount)
{
}

void SceneNode::flatten()
{
  for( NodeList::const_iterator it = m_children.begin();
       it != m_children.end();
       it++ )
  {
    (*it)->m_trans = m_trans * (*it)->m_trans;
    (*it)->flatten();
  }
}

void SceneNode::print(int depth) const
{
  std::cerr << std::string(depth+depth, ' ') << m_name << std::endl;
  for( NodeList::const_iterator it = m_children.begin();
       it != m_children.end();
       it++ )
  {
    (*it)->print(depth+1);
  }
}

// GeometryNode

GeometryNode::GeometryNode(const std::string& name, Primitive* primitive)
  : SceneNode(name)
  , m_primitive(primitive)
  , m_material(NULL)
{
}

GeometryNode::GeometryNode(const aiMesh *mesh)
  : SceneNode(mesh->mName.C_Str())
  , m_primitive(new Mesh(mesh))
  , m_material(NULL)
{
}

GeometryNode::GeometryNode(const GeometryNode &orig)
  : SceneNode(orig)
  , m_primitive(orig.m_primitive)
  , m_material(orig.m_material)
{
}

GeometryNode::~GeometryNode()
{ 
  // To make this node disown its members, set them to NULL.
  // This way we may transfer ownership elsewhere.
  if (m_material)
    delete m_material;

  if (m_primitive)
    delete m_primitive; 
}

void GeometryNode::print(int depth) const
{
  if (m_primitive->is_mesh())
  {
    std::cerr << *(Mesh *)m_primitive << std::cerr;
  }
  SceneNode::print(depth);
}
