#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include <boost/shared_ptr.hpp>
#include "eigen.h"

class Primitive
{
public:
  Primitive()
    : m_bbox(Vector3f(-1.0f,-1.0f,-1.0f), Vector3f(1.0f,1.0f,1.0f))
  { }
  virtual ~Primitive() { qDebug() << "destroying primitive:" << this; }
  virtual inline bool is_mesh()       const { return false; }
  virtual inline bool is_pointcloud() const { return false; }
  virtual inline bool is_dynamic()    const { return false; }
  virtual inline void transform_in_place(const AffineCompact3f &trans) = 0;
  AlignedBox3f &get_bbox() { return m_bbox; }
  virtual inline AlignedBox3f &compute_bbox() = 0;
  inline void cube_bbox() { m_bbox = AlignedBox3f(Vector3f(-1.f,-1.f,-1.f), Vector3f(1.f,1.f,1.f)); }
protected:
  AlignedBox3f m_bbox;
};

typedef boost::shared_ptr<Primitive> PrimitivePtr;

#endif // PRIMITIVE_H
