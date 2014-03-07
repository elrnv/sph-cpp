#ifndef MATERIAL_H
#define MATERIAL_H

#include <assimp/material.h>
#include "eigen.h"

// Materials

class Material {
public:
  Material(const aiMaterial &mat);
  Material(const Vector3f &ka, const Vector3f &kd, const Vector3f &ks,
                float shininess,
                float reflectivity);
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

  float m_opacity;
  float m_shininess;
  float m_reflectivity;
};

// define a default material
const Material DEFAULT_MATERIAL =
  Material(Vector3f(0.0f,0.0f,0.0f),
           Vector3f(0.6f,0.6f,0.6f),
           Vector3f(0.2f,0.2f,0.0f),
           25.0f, 0.0f);

#endif
