#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include "eigen.h"

class Primitive {
public:
  Primitive() { }
  virtual ~Primitive() { }
  virtual bool is_mesh() const { return false; }
};

class Sphere : public Primitive 
{
  // The default sphere is centered at the origin with radius 1

public:
  Sphere() { }
};

class Cylinder: public Primitive
{
  // The default cylinder is hierarchical and is positioned at the origin with
  // radius 1 and height 2
  Vector3d m_normal_t; // normal of the top cap (bottom cap is the negative)

public:
  Cylinder()
    : m_normal_t(0,1,0)
  { }
  
};

class Cone: public Primitive
{
  // The default cone is hierarchical with its base positioned at the origin
  // and the tip distance 1 above the origin
  Vector3d m_normal_b;

public:
  Cone()
    : m_normal_b(0,-1,0) 
  { }
};

#endif
