#include "material.h"

// Material

Material::~Material()
{
}

PhongMaterial::PhongMaterial(const Color& kd, const Color& ks,
                             double shininess, double reflectivity)
  : m_kd(kd), m_ks(ks), m_shininess(shininess), m_reflectivity(reflectivity)
{
}

PhongMaterial::~PhongMaterial()
{
}

