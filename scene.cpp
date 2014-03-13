#include <iostream>
#include <limits>
#include "eigen.h"
#include "scene.h"
#include "mesh.h"
#include "pointcloud.h"

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

  Matrix4f mat;
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
  for( const SceneNode *node : orig.m_children )
  {
    SceneNode *new_node = node->clone(); // dynamic allocation with new
    m_children.push_back(new_node);
  }
}

SceneNode::~SceneNode()
{
}

void SceneNode::rotate(float angle, const Vector3f &axis)
{
  m_trans.rotate(AngleAxisf(angle*RADIAN, axis));
}

void SceneNode::scale(float amount)
{
  m_trans.scale(amount);
}

void SceneNode::scale(const Vector3f& amount)
{
  m_trans.scale(amount);
}

void SceneNode::translate(const Vector3f& amount)
{
  m_trans.translate(amount);
}

void SceneNode::flatten()
{
  for( SceneNode *node : m_children )
  {
    node->m_trans = m_trans * node->m_trans;
    node->flatten();
  }
  m_trans.setIdentity();
}

void SceneNode::print(int depth) const
{
  std::cerr << std::string(depth+depth, ' ') << m_name << std::endl;
  for( const SceneNode *node : m_children )
    node->print(depth+1);
}

unsigned int SceneNode::num_primitives() const
{
  unsigned int count = 0;
  for( const SceneNode *node : m_children )
    count += node->num_primitives();

  return count + (is_geometry() ? 1 : 0);
}

AlignedBox3f &SceneNode::compute_bbox()
{
  m_bbox.setEmpty();
  for( SceneNode *node : m_children )
    m_bbox.extend(node->compute_bbox());
  return m_bbox;
}
void SceneNode::cube_bbox()
{ 
  m_bbox = AlignedBox3f(Vector3f(-1.f,-1.f,-1.f), Vector3f(1.f,1.f,1.f)); 
  for( SceneNode *node : m_children )
    node->cube_bbox();
}

void SceneNode::normalize_model()
{
  normalize_model(Vector3f(0.0f, 0.0f, 0.0f),Vector3f(0.0f, 0.0f, 0.0f));
}


void SceneNode::normalize_model(
    const Vector2f &ext_x,
    const Vector2f &ext_y,
    const Vector2f &ext_z)
{
  normalize_model(Vector3f(ext_x[0], ext_y[0], ext_z[0]),
                  Vector3f(ext_x[1], ext_y[1], ext_z[1]));
}

void SceneNode::normalize_model(const Vector3f &ext_blf, const Vector3f &ext_trc)
{
  compute_bbox();
  Vector3f blf( m_bbox.corner(Eigen::AlignedBox3f::BottomLeftFloor) );
  Vector3f trc( m_bbox.corner(Eigen::AlignedBox3f::TopRightCeil) );
  m_bbox.extend( blf - ext_blf );
  m_bbox.extend( trc + ext_trc );

  Vector3f sizevec = m_bbox.sizes();
  m_trans.scale(2.0f/sizevec.maxCoeff());
  Vector3f box_center = m_bbox.center();
  m_trans.translate(-box_center);
}




// GeometryNode

GeometryNode::GeometryNode(const std::string& name, Primitive* primitive)
  : SceneNode(name)
  , m_primitive(primitive)
  , m_material(&DEFAULT_MATERIAL)
{
}

GeometryNode::GeometryNode(const aiMesh *mesh, const aiMaterial *mat)
  : SceneNode(mesh->mName.C_Str())
  , m_primitive(NULL)
  , m_material(&DEFAULT_MATERIAL)
{
  if ( mesh->mPrimitiveTypes & aiPrimitiveType_POINT )
  {
    m_primitive = new PointCloud(mesh); // interpret as point cloud
  }
  else if ( mesh->mPrimitiveTypes & aiPrimitiveType_TRIANGLE)
  {
    m_primitive = new Mesh(mesh); // interpret as triangular mesh
  }
  // otherwise m_primitive remains NULL

  if (!mat)
    return;

  m_material = new Material(*mat);
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

void GeometryNode::flatten()
{
  if (m_primitive)
    m_primitive->transform_in_place(m_trans);

  SceneNode::flatten();
}

void GeometryNode::print(int depth) const
{
  if (!m_primitive)
  {
    std::cerr << "NULL Primitive." << std::endl;
  }
  else if (m_primitive->is_mesh())
  {
    std::cerr << *(Mesh *)m_primitive << std::endl;
  }
  SceneNode::print(depth);
}

AlignedBox3f &GeometryNode::compute_bbox()
{
  m_bbox.setEmpty();
  if (m_primitive)
    m_bbox.extend(m_primitive->compute_bbox());
  for( SceneNode *node : m_children )
    m_bbox.extend(node->compute_bbox());
  return m_bbox;
}

void GeometryNode::cube_bbox()
{
  SceneNode::cube_bbox();
  if (m_primitive)
    m_primitive->cube_bbox();
}
