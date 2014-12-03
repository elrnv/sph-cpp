#ifndef MATERIAL_H
#define MATERIAL_H

#include <assimp/material.h>
#include <boost/shared_ptr.hpp>
#include "eigen.h"

// Materials

class Material {
public:
  Material(float opacity = 0.3f); // default material
  Material(const Material &mat);
  Material(const aiMaterial &mat);
  Material(const Vector3f &ka, const Vector3f &kd, const Vector3f &ks,
                float shininess, float reflectivity, float opacity);
  virtual ~Material();

  const Vector3f &ka() const { return m_ka; }
  const Vector3f &kd() const { return m_kd; }
  const Vector3f &ks() const { return m_ks; }
  float get_opacity() const { return m_opacity; }
  float get_shininess() const { return m_shininess; }
  float get_reflectivity() const { return m_reflectivity; }

private:
  Vector3f m_ka; // ambient
  Vector3f m_kd; // diffuse
  Vector3f m_ks; // specular

  float m_shininess;
  float m_reflectivity;
  float m_opacity;
};

typedef boost::shared_ptr< const Material > MaterialConstPtr;

#endif
