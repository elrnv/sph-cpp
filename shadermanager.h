#ifndef SHADERMANAGER_H
#define SHADERMANAGER_H
#include <QtGui/QOpenGLShaderProgram>

class OpenGLWindow;

class ShaderManager : public QObject
{
public:
  ShaderManager(OpenGLWindow *parent);
  ~ShaderManager();

  enum ShaderType
  {
    NORMALS,
    FLAT,
    PHONG,
    PARTICLE,
    ADDITIVE_PARTICLE
  };

  QOpenGLShaderProgram *get_normals_shader() { return &m_normals_shader; }
  QOpenGLShaderProgram *get_flat_shader() { return &m_flat_shader; }
  QOpenGLShaderProgram *get_phong_shader()     { return &m_phong_shader; }
  QOpenGLShaderProgram *get_particle_shader()  { return &m_particle_shader; }
  QOpenGLShaderProgram *get_additive_particle_shader()  { return &m_add_particle_shader; }

  void init();

private:
  QOpenGLShaderProgram m_normals_shader;
  QOpenGLShaderProgram m_flat_shader;
  QOpenGLShaderProgram m_phong_shader;
  QOpenGLShaderProgram m_particle_shader;
  QOpenGLShaderProgram m_add_particle_shader;
};

#endif // SHADERMANAGER_H
