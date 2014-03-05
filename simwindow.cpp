#include <QtGui/QScreen>
#include <QtGui/QKeyEvent>
#include <string>
#include "scene.h"
#include "util.h"
#include "mesh.h"
#include "simwindow.h"

SimWindow::SimWindow()
  : m_viewmode(ShaderManager::PHONG)
  , m_shaderman(this)
{

}

void SimWindow::init()
{
  m_shaderman.init();

	m_global_uniform.create();
	m_global_uniform.setUsagePattern( UniformBuffer::StreamDraw );
	m_global_uniform.bind();
  m_global_uniform.allocate( sizeof(m_ubo) );
  m_global_uniform.bindToIndex();

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
  update_viewmode(m_viewmode);
  load_model(1);
}

void SimWindow::load_model(int i)
{
  std::string filename;

  // optionally extend centered container box on load
  Vector2f ext_x(0.0f, 0.0f);
  Vector2f ext_y(0.0f, 0.0f);
  Vector2f ext_z(0.0f, 0.0f);

  // optionally rotate (useful to change coordinates)
  double angle_x = 0.0f;
  double angle_y = 0.0f;

  // optionally uniformly scale model
  float scale = 1.0f;


  switch (i)
  {
    case 0:
      filename = "test.obj";
      break;
    case 1:
      filename = "cube.obj";
      break;
    case 2:
      filename = "cow.obj";
      angle_y = 180.0f;
      break;
    case 3:
      filename = "bun_zipper.ply";
      break;
    case 4:
      filename = "bunnyData.obj";
      break;
    case 5:
      filename = "pcm-house4.obj";
      angle_x = -90.0f;
      break;
    case 6:
      filename = "cowpts.obj";
      angle_y = 180.0f;
      break;
    case 7:
      filename = "sphere.obj";
      ext_x = Vector2f(2.0f, 2.0f);
      ext_y = Vector2f(4.0f, 0.0f);
      ext_z = Vector2f(2.0f, 2.0f);
      scale = 0.7f;
      break;
  }
  SceneNode *scene = Util::loadScene(filename);
  if (!scene)
    return;

  scene->scale(scale);
  scene->rotate(angle_x, Vector3f::UnitX());
  scene->rotate(angle_y, Vector3f::UnitY());
  scene->normalize_model(ext_x, ext_y, ext_z);
  m_ubo.modelmtx = Affine3d(scene->get_trans()).matrix().cast<float>();

#ifndef QT_NO_DEBUG
  scene->print();
#endif
  
  m_glprims.clear();
  m_glprims.reserve( scene->num_primitives() );
  Util::loadGLData( scene, m_glprims, m_global_uniform, m_shaderman );
  delete scene;
}

void SimWindow::update_viewmode(ViewMode vm)
{
  m_viewmode = vm;

  glEnable(GL_BLEND);

  if (vm == ViewMode::PARTICLE) 
  {
    glDisable(GL_DEPTH_TEST);
	  glDisable(GL_CULL_FACE);
    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_PROGRAM_POINT_SIZE);
  }
  else
  {
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_PROGRAM_POINT_SIZE);
  }
}

void SimWindow::make_dynamic()
{
  for ( const GLPrimitivePtr &prim_ptr : m_glprims )
  {
    GLPrimitive *glprim = prim_ptr.get();
    if (!glprim->is_pointcloud())
      continue;
  
    GLPointCloud *glpc = static_cast<GLPointCloud *>(glprim);
    glpc->make_dynamic(/*particle mass = */1.0f);
  }
  set_animating(true);
}

void SimWindow::render()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

  for ( const GLPrimitivePtr &prim_ptr : m_glprims )
  {
    GLPrimitive *glprim = prim_ptr.get();

    Affine3f vtrans(get_view_trans()); // view transformation
    m_ubo.vpmtx = (get_proj_trans() * vtrans).matrix();
    m_ubo.mvpmtx = m_ubo.vpmtx * m_ubo.modelmtx;
    m_ubo.vinvmtx = vtrans.inverse().matrix();
    m_ubo.eyepos = m_ubo.vinvmtx * Vector4f(0.0f, 0.0f, 0.0f, 1.0f);

    m_ubo.normalmtx.block(0,0,3,3) = m_ubo.modelmtx.block(0,0,3,3).inverse().transpose();

    glprim->update_shader(m_viewmode);

    glprim->get_program()->bind();
    m_global_uniform.bind();

    int offset = 0;
    offset = m_global_uniform.write(offset, m_ubo.mvpmtx.data(), sizeof(Matrix4f));
    offset = m_global_uniform.write(offset, m_ubo.vpmtx.data(), sizeof(Matrix4f));
    offset = m_global_uniform.write(offset, m_ubo.modelmtx.data(), sizeof(Matrix4f));
    offset = m_global_uniform.write(offset, m_ubo.normalmtx.data(), sizeof(Matrix4f));
    offset = m_global_uniform.write(offset, m_ubo.vinvmtx.data(), sizeof(Matrix4f));
    offset = m_global_uniform.write(offset, m_ubo.eyepos.data(), sizeof(Vector4f));
    
    Vector4f l1 = m_ubo.vinvmtx * Vector4f(0.0, 0.0, 10.0, 1.0);
    glprim->get_program()->setUniformValue("lights[0].pos", QVector4D(l1[0], l1[1], l1[2], l1[3]));
    glprim->get_program()->setUniformValue("lights[0].col", QVector4D(0.8, 0.9, 0.9, 1.0));
    glprim->get_program()->setUniformValue("lights[1].pos", QVector4D(0.0, 0.0, 0.0, 0.0));
    glprim->get_program()->setUniformValue("lights[1].col", QVector4D(0.0, 0.0, 0.0, 0.0));
    glprim->get_program()->setUniformValue("ambientColor", QVector4D(0.01, 0.02, 0.01, 1.0));
    glprim->get_program()->setUniformValue("diff", QVector4D(0.2, 0.5, 1.0, 1.0));
    glprim->get_program()->setUniformValue("spec", QVector4D(0.5, 1.0, 1.0, 1.0));
    glprim->get_program()->setUniformValue("specpow", 25.0f);

    glprim->get_vao().bind();

    if (m_viewmode == ShaderManager::PARTICLE)
    {
      glprim->get_program()->setUniformValue("pt_scale", float(14*window_dim()[1]*m_near));
      glprim->get_program()->setUniformValue("pt_radius", 0.03f);
      glDrawArrays(GL_POINTS, 0, glprim->get_num_vertices());
    }
    else
    {
      glDrawElements(GL_TRIANGLES, glprim->get_num_indices(), GL_UNSIGNED_INT, 0);
    }

    glprim->get_vao().release();
    m_global_uniform.release();
    glprim->get_program()->release();
  }
}

void SimWindow::keyPressEvent(QKeyEvent *event)
{
  int key = event->key();
  switch (key)
  {
    case Qt::Key_D:
      make_dynamic();
      break;
    case Qt::Key_W:
      update_viewmode(ViewMode::WIREFRAME);
      break;
    case Qt::Key_S:
      update_viewmode(ViewMode::PHONG);
      break;
    case Qt::Key_P:
      update_viewmode(ViewMode::PARTICLE);
      break;
    case Qt::Key_1: load_model(1); break;
    case Qt::Key_2: load_model(2); break;
    case Qt::Key_3: load_model(3); break;
    case Qt::Key_4: load_model(4); break;
    case Qt::Key_5: load_model(5); break;
    case Qt::Key_6: load_model(6); break;
    case Qt::Key_7: load_model(7); break;
    case Qt::Key_8: load_model(8); break;
    case Qt::Key_9: load_model(9); break;
    case Qt::Key_0: load_model(0); break;
    default:
      OpenGLWindow::keyPressEvent(event);
      return;
  }

  renderLater(); // queue rendering
}

