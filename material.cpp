#include "material.h"

// Material

// define the default material
Material::Material()
 : m_ka(0.0f,0.0f,0.0f)
 , m_kd(0.6f,0.6f,0.6f)
 , m_ks(0.2f,0.2f,0.0f)
 , m_shininess(25.0f)
 , m_reflectivity(0.0f)
 , m_opacity(0.3f)
{}

// define the copy constructor
Material::Material(const Material &orig)
 : m_ka(orig.m_ka)
 , m_kd(orig.m_kd)
 , m_ks(orig.m_ks)
 , m_shininess(orig.m_shininess)
 , m_reflectivity(orig.m_reflectivity)
 , m_opacity(orig.m_opacity)
{}

// material imported from assimp
Material::Material(const aiMaterial &mat)
{
  aiColor3D ka, kd, ks;
  if (AI_SUCCESS != mat.Get<aiColor3D>(AI_MATKEY_COLOR_AMBIENT, ka))
    ka = aiColor3D(0.0f, 0.0f, 0.0f);
  if (AI_SUCCESS != mat.Get<aiColor3D>(AI_MATKEY_COLOR_DIFFUSE, kd))
    kd = aiColor3D(0.6f, 0.6f, 0.6f);
  if (AI_SUCCESS != mat.Get<aiColor3D>(AI_MATKEY_COLOR_SPECULAR, ks))
    ks = aiColor3D(0.0f, 0.0f, 0.0f);
  if (AI_SUCCESS != mat.Get<float>(AI_MATKEY_SHININESS, m_shininess))
    m_shininess = 0.f;

  if (AI_SUCCESS != mat.Get<float>(AI_MATKEY_OPACITY, m_opacity))
    m_opacity = 0.3f;

  m_ka = Vector3f(ka.r, ka.g, ka.b);
  m_kd = Vector3f(kd.r, kd.g, kd.b);
  m_ks = Vector3f(ks.r, ks.g, ks.b);
  m_reflectivity = 0.0f;
}

// custom material
Material::Material(const Vector3f& ka, const Vector3f& kd, const Vector3f& ks,
    float shininess, float reflectivity, float opacity)
  : m_ka(ka), m_kd(kd), m_ks(ks)
  , m_shininess(shininess), m_reflectivity(reflectivity), m_opacity(opacity)
{
}

Material::~Material()
{ }

