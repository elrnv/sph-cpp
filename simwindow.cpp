#include <QtGui/QScreen>
#include "scene.h"
#include "util.h"
#include "mesh.h"
#include "simwindow.h"

SimWindow::SimWindow()
  : m_num_vertices(0)
  , m_num_indices(0)
	, m_posBuffer(QOpenGLBuffer::VertexBuffer)
	, m_nmlBuffer(QOpenGLBuffer::VertexBuffer)
	, m_idxBuffer(QOpenGLBuffer::IndexBuffer)
	, m_program(0)
	, m_frame(0)
{

}

void SimWindow::init()
{
	m_program = new QOpenGLShaderProgram(this);
#if 0
	m_program->addShaderFromSourceFile(QOpenGLShader::Vertex, ":/normals.vert");
	m_program->addShaderFromSourceFile(QOpenGLShader::Geometry, ":/normals.geom");
	m_program->addShaderFromSourceFile(QOpenGLShader::Fragment, ":/plain.frag");
#else
	m_program->addShaderFromSourceFile(QOpenGLShader::Vertex, ":/phong.vert");
	m_program->addShaderFromSourceFile(QOpenGLShader::Fragment, ":/phong.frag");
#endif
  qDebug() << "Shader Program Log:" << m_program->log();
	m_program->link();

  SceneNode *scene = Util::loadScene("cube.obj");
  scene->print();

  Util::extractGLSceneDataSize<GLuint>(scene, m_num_vertices, m_num_indices);

  GLfloat vertices[m_num_vertices];
	GLfloat normals[m_num_vertices];
	GLubyte indices[m_num_indices];
  qDebug() << m_num_indices;
  Util::extractGLSceneData<GLfloat, GLubyte>(scene, vertices, normals, indices);

  for (int i = 0; i < m_num_vertices; ++i)
    qDebug() << normals[i];

	m_vao1 = new QOpenGLVertexArrayObject( this );
	m_vao1->create();
	m_vao1->bind();

	m_posBuffer.create();
	m_posBuffer.setUsagePattern( QOpenGLBuffer::StaticDraw );
	m_posBuffer.bind();
	m_posBuffer.allocate( vertices, sizeof( vertices ) );
	m_program->enableAttributeArray( "pos" );
	m_program->setAttributeBuffer( "pos", GL_FLOAT, 0, 3 );

	m_nmlBuffer.create();
	m_nmlBuffer.setUsagePattern( QOpenGLBuffer::StaticDraw );
	m_nmlBuffer.bind();
	m_nmlBuffer.allocate( normals, sizeof( normals ) );
	m_program->enableAttributeArray( "nml" );
	m_program->setAttributeBuffer( "nml", GL_FLOAT, 0, 3 );

	m_idxBuffer.create();
	m_idxBuffer.setUsagePattern( QOpenGLBuffer::StaticDraw );
	m_idxBuffer.bind();
	m_idxBuffer.allocate( indices, sizeof( indices ) );

	m_mvpMtxUniform = m_program->uniformLocation("mvpMtx");
	m_mvInvMtxUniform = m_program->uniformLocation("mvInvMtx");
	m_nmlMtxUniform = m_program->uniformLocation("nmlMtx");
	m_mvMtxUniform = m_program->uniformLocation("mvMtx");

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
}

void SimWindow::render()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	m_program->bind();

	m_vao1->bind();

	AffineCompact3f MV = get_view_mtx(); // view * model transformation
	AffineCompact3f MVinv = MV.inverse(Eigen::Affine); // view * model transformation

	m_program->setUniformValue(m_mvInvMtxUniform, MVinv.toQMatrix4x4());
	m_program->setUniformValue(m_mvpMtxUniform, (get_proj_mtx() * MV).toQMatrix4x4());
	m_program->setUniformValue(m_nmlMtxUniform, QMatrix3x3(MVinv.linear().data())); // Eigen matrices are ColMajor
	m_program->setUniformValue(m_mvMtxUniform, MV.toQMatrix4x4()); // Eigen matrices are ColMajor

  m_program->setUniformValue("lights[0].pos", QVector4D(5.0, 5.0, 5.0, 0.0));
  m_program->setUniformValue("lights[0].col", QVector4D(1.0, 0.5, 0.5, 0.5));
  m_program->setUniformValue("lights[0].camera", true);
  m_program->setUniformValue("ambientMat", QVector4D(0.01, 0.02, 0.01, 1.0));
  m_program->setUniformValue("diffuseMat", QVector4D(0.5, 1.0, 0.5, 1.0));
  m_program->setUniformValue("specMat", QVector4D(0.5, 1.0, 0.5, 1.0));
  m_program->setUniformValue("specPow", 25.0f);

	//m_program->setUniformValue(m_colUniform, QVector3D(1.0f, 0.0f, 0.0f));

	glDrawElements(GL_TRIANGLES, m_num_indices, GL_UNSIGNED_BYTE, 0);

	m_program->release();

	++m_frame;
}
