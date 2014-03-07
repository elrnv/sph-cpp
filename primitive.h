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
  virtual inline void transform(const AffineCompact3f &trans) = 0;
  AlignedBox3f &get_bbox() { return m_bbox; }
protected:
  AlignedBox3f m_bbox;
};

#endif // PRIMITIVE_H
