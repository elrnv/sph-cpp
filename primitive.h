#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include "eigen.h"

class Primitive
{
public:
  Primitive()
    : m_bbox(Vector3f(-1.0f,-1.0f,-1.0f), Vector3f(1.0f,1.0f,1.0f))
  { }
  virtual ~Primitive() { }
  virtual inline bool is_mesh()       const { return false; }
  virtual inline bool is_pointcloud() const { return false; }
  virtual inline bool is_dynamic()    const { return false; }
  AlignedBox3f &get_bbox() { return m_bbox; }
protected:
  AlignedBox3f m_bbox;
};

class Sphere : public Primitive 
{
public:
  Sphere()
  { }
};

class Cylinder : public Primitive
{
public:
  Cylinder()
    : m_normal_t(0,1,0)
  { }
  
private:
  // The default cylinder is hierarchical and is positioned at the origin with
  // radius 1 and height 2
  Vector3d m_normal_t; // normal of the top cap (bottom cap is the negative)

};

class Cone : public Primitive
{
public:
  Cone()
    : m_normal_b(0,-1,0) 
  { }

private:
  // The default cone is hierarchical with its base positioned 1 below the origin
  // and the tip distance 1 above the origin, that is pointing up
  Vector3d m_normal_b;
};

#endif // PRIMITIVE_H
