#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include <boost/shared_ptr.hpp>
#include "materialmanager.h"
#include "eigen.h"

class Primitive
{
public:
  Primitive(Index matidx)
    : m_material_idx(matidx)
    , m_bbox(Vector3f(-1.0f,-1.0f,-1.0f), Vector3f(1.0f,1.0f,1.0f))
  { }

  virtual ~Primitive() { }

  inline Index get_material_idx() const       { return m_material_idx; }
  inline void  set_material_idx(Index matidx) { m_material_idx = matidx; }
  inline Vector3f get_color(MaterialManager &matman) const 
  { return matman[m_material_idx].kd(); }

  virtual bool is_mesh()       const = 0;
  virtual bool is_pointcloud() const = 0;
  virtual inline void transform_in_place(const AffineCompact3f &trans) = 0;
  AlignedBox3f &get_bbox() { return m_bbox; }
  virtual inline AlignedBox3f &compute_bbox() = 0;
  inline void cube_bbox() 
  { 
    m_bbox = AlignedBox3f(Vector3f(-1.f,-1.f,-1.f), Vector3f(1.f,1.f,1.f)); 
  }
protected:
  Index        m_material_idx; // index of the material being used
  AlignedBox3f m_bbox;
};

typedef boost::shared_ptr<Primitive> PrimitivePtr;

#endif // PRIMITIVE_H
