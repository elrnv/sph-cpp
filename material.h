#ifndef MATERIAL_H
#define MATERIAL_H

#include "eigen.h"

// Materials

class Material {
public:
  virtual ~Material();

protected:
  Material()
  {
  }
};

class PhongMaterial : public Material {
public:
  PhongMaterial(const Color& kd, const Color& ks,
                double shininess,
                double reflectivity);
  virtual ~PhongMaterial();

  const Color &get_kd() const { return m_kd; }
  const Color &get_ks() const { return m_ks; }
  double get_shininess() const { return m_shininess; }
  double get_reflectivity() const { return m_reflectivity; }

private:
  Color m_kd;
  Color m_ks;

  double m_shininess;
  double m_reflectivity;
};

// define a default material
const PhongMaterial DEFAULT_MATERIAL =
  PhongMaterial(Color(0.6f,0.6f,0.6f,1.0f), Color(0.0f,0.0f,0.0f,1.0f), 25.0f, 0.0f);

#endif
