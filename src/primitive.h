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

  virtual bool is_mesh()       const { return false; }
  virtual bool is_pointcloud() const { return false; }
  virtual inline void transform_in_place(const Affine3f &trans) = 0;
  virtual inline AlignedBox3f &compute_bbox() = 0;

  inline AlignedBox3f &get_bbox() { return m_bbox; }
  inline void set_bbox(const AlignedBox3f &bbox) { m_bbox = bbox; }

protected:
  Index        m_material_idx; // index of the material being used
  AlignedBox3f m_bbox;
};

typedef boost::shared_ptr<Primitive> PrimitivePtr;

#endif // PRIMITIVE_H
