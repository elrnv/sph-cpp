#include "openglwindow.h"
#include "shadermanager.h"

ShaderManager::ShaderManager(OpenGLWindow *parent)
  : QObject(parent)
  , m_wireframe_shader(this)
  , m_phong_shader(this)
  , m_particle_shader(this)
{
}
ShaderManager::~ShaderManager()
{

}

void ShaderManager::init()
{
  m_wireframe_shader.addShaderFromSourceFile(QOpenGLShader::Vertex,   ":/normals.vert");
  m_wireframe_shader.addShaderFromSourceFile(QOpenGLShader::Geometry, ":/normals.geom");
  m_wireframe_shader.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/plain.frag");
  m_wireframe_shader.link();
  qDebug() << "Wireframe Shader LOG:" << m_wireframe_shader.log();

  m_phong_shader.addShaderFromSourceFile(QOpenGLShader::Vertex,   ":/phong.vert");
  m_phong_shader.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/phong.frag");
  m_phong_shader.link();
  qDebug() << "Phong Shader LOG:" << m_phong_shader.log();

  m_particle_shader.addShaderFromSourceFile(QOpenGLShader::Vertex, ":/particle.vert");
  m_particle_shader.addShaderFromSourceFile(QOpenGLShader::Geometry, ":/particle.geom");
  m_particle_shader.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/particle.frag");
  m_particle_shader.link();
  qDebug() << "Particle Shader LOG:" << m_particle_shader.log();
}
