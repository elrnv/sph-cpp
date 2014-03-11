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
    WIREFRAME,
    PHONG,
    PARTICLE
  };

  QOpenGLShaderProgram *get_wireframe_shader() { return &m_wireframe_shader; }
  QOpenGLShaderProgram *get_phong_shader()     { return &m_phong_shader; }
  QOpenGLShaderProgram *get_particle_shader()  { return &m_particle_shader; }

  void init();

private:
  QOpenGLShaderProgram m_wireframe_shader;
  QOpenGLShaderProgram m_phong_shader;
  QOpenGLShaderProgram m_particle_shader;
};

#endif // SHADERMANAGER_H
