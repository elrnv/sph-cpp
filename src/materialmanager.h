#ifndef MATERIALMANAGER_H
#define MATERIALMANAGER_H

#include <iostream>
#include "material.h"
#include "types.h"

// MaterialManager keeps track of materials used throughout the scene

class MaterialManager
{
public:
  MaterialManager()
    : m_materials(1, Material())
  {
    // also add a default fully transparent material
    m_materials.push_back(Material(0.0f));
  }
  ~MaterialManager() { }

  Index get_transparent_material_idx() { return 1; }
  Index get_default_material_idx() { return 0; }

  /// Manager interface
  Index add_material(const aiMaterial *mat)
  {
    if (!mat)
      return 0; // default material

    m_materials.push_back( Material(*mat) );
    return m_materials.size() - 1;
  }

  typedef std::vector< Material > MaterialVec;

  inline const Material &get_material(Index matidx) const
  { 
    return m_materials[matidx]; 
  }

  inline const Material &operator[] (Index matidx) const
  { 
    return m_materials[matidx]; 
  }

  inline Material &operator[] (Index matidx)
  { 
    return m_materials[matidx]; 
  }

private:
  MaterialVec m_materials;
};

#endif // MATERIALMANAGER_H
