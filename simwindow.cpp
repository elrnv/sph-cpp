#include <QtGui/QScreen>
#include <QtGui/QKeyEvent>
#include <string>
#include "scene.h"
#include "util.h"
#include "mesh.h"
#include "simwindow.h"
#include "dynamics.h"

void SimWindow::toggle_shortcuts()
{
  m_show_shortcuts = !m_show_shortcuts;
  if (!m_show_shortcuts)
  {
    glclear_bl();
    return;
  }

  glprintf_bl("Shortcuts:\n");
  glprintf_blc(BLUE, "  H: show/hide shortcuts\n");
  glprintf_bl("  View Modes: \n");
  glprintf_blc(GREEN, "    W: wireframe\n");
  glprintf_blc(RED,   "    S: phong\n");
  glprintf_blc(BLUE, "    P: particles\n");
  glprintf_bl("  Models: \n");
  glprintf_blc(GREEN, "    1: cube\n");
  glprintf_blc(GREEN, "    2: cow\n");
  glprintf_blc(GREEN, "    3: bunny\n");
  glprintf_blc(BLUE,  "    4: sparse point cloud\n");
  glprintf_blc(BLUE,  "    5: normal point cloud\n");
  glprintf_blc(BLUE,  "    6: dense point cloud\n");
  glprintf_bl("  Dynamics: \n");
  glprintf_blc(RED,  "    D: enable\n");
  glprintf_blc(RED,  "    C: disable\n");
}

SimWindow::SimWindow()
  : m_show_shortcuts(false) // immediately toggled below
  , m_viewmode(ShaderManager::PARTICLE)
  , m_change_prog(true)
  , m_shaderman(this)
{
  toggle_shortcuts();
}

SimWindow::~SimWindow()
{
  clear_threads();
}

void SimWindow::clear_threads()
{
  for ( const GLPrimitivePtr &prim_ptr : m_glprims )
  {
    GLPrimitive *glprim = prim_ptr.get();
    if (!glprim->is_pointcloud())
      continue;

    GLPointCloud *glpc = static_cast<GLPointCloud *>(glprim);
    PointCloud *pc = glpc->get_pointcloud();

    if (!pc->is_dynamic())
      continue;

    DynamicPointCloud *dpc = static_cast<DynamicPointCloud *>(pc);
    dpc->request_stop(); // stop threads
  }

  for ( auto & thread : m_sim_threads )
    thread.join();

  m_sim_threads.clear();
}

void SimWindow::init()
{
  OpenGLWindow::init();
  m_shaderman.init();

	m_ubo.create();
	m_ubo.setUsagePattern( UniformBuffer::StreamDraw );
	m_ubo.bind();
  m_ubo.allocate( sizeof(m_udata) );
  m_ubo.bindToIndex();

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
  load_model(4);
}

void SimWindow::load_model(int i)
{
  clear_threads();
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
      filename = "sparsesphere.obj";
      ext_x = Vector2f(1.0f, 1.0f);
      ext_y = Vector2f(2.0f, 0.0f);
      ext_z = Vector2f(1.0f, 1.0f);
      angle_x = 40.0f;
      angle_y = 40.0f;
      break;
    case 5:
      filename = "sphere.obj";
      ext_x = Vector2f(1.0f, 1.0f);
      ext_y = Vector2f(2.0f, 0.0f);
      ext_z = Vector2f(1.0f, 1.0f);
      angle_x = 40.0f;
      angle_y = 40.0f;
      break;
    case 6:
      filename = "densesphere.obj";
      ext_x = Vector2f(1.0f, 1.0f);
      ext_y = Vector2f(2.0f, 0.0f);
      ext_z = Vector2f(1.0f, 1.0f);
      angle_x = 40.0f;
      angle_y = 40.0f;
      break;
    default:
      glprintf_trc(RED, "model not found\n");
      return;
  }
  SceneNode *scene = Util::loadScene(filename);
  if (!scene)
    return;

  scene->scale(scale);
  scene->rotate(angle_x, Vector3f::UnitX());
  scene->rotate(angle_y, Vector3f::UnitY());
  scene->normalize_model(ext_x, ext_y, ext_z);
  scene->flatten();
  scene->cube_bbox();
  m_udata.modelmtx.setIdentity();

#ifndef QT_NO_DEBUG
  scene->print();
#endif
  
  m_glprims.clear();
  m_glprims.reserve( scene->num_primitives() );
  Util::loadGLData( scene, m_glprims, m_ubo, m_shaderman );
  delete scene;
  change_viewmode(m_viewmode);
}

void SimWindow::change_viewmode(ViewMode vm)
{
  m_viewmode = vm;
  m_change_prog = true;
  reset_viewmode();
}

void SimWindow::reset_viewmode()
{
  glEnable(GL_BLEND);

  if (m_viewmode == ViewMode::PARTICLE) 
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

void SimWindow::make_static()
{
  clear_threads();
  set_animating(false);
}

void SimWindow::make_dynamic()
{
  clear_threads();
  for ( const GLPrimitivePtr &prim_ptr : m_glprims )
  {
    GLPrimitive *glprim = prim_ptr.get();
    if (!glprim->is_pointcloud())
      continue;
  
    // use right text buffer for dynamic properties
    glclear_tr(); // so clear it

    GLPointCloud *glpc = static_cast<GLPointCloud *>(glprim);
    DynamicPointCloud *dpc = glpc->make_dynamic(
        /*density = */1000.0f,
        /*viscosity = */1000.0f,
        /*surface tension coefficient = */0.0728f);

    // run simulation
    m_sim_threads.push_back(std::thread(&DynamicPointCloud::run, dpc));
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
    m_udata.vpmtx = (get_proj_trans() * vtrans).matrix();
    m_udata.mvpmtx = m_udata.vpmtx * m_udata.modelmtx;
    m_udata.vinvmtx = vtrans.inverse().matrix();
    m_udata.eyepos = m_udata.vinvmtx * Vector4f(0.0f, 0.0f, 0.0f, 1.0f);

    m_udata.normalmtx.block(0,0,3,3) = m_udata.modelmtx.block(0,0,3,3).inverse().transpose();

    m_udata.ambient = glprim->get_ambient();
    m_udata.diffuse = glprim->get_diffuse();
    m_udata.specular = glprim->get_specular();
    m_udata.options[0] = glprim->get_specpow();
    m_udata.options[1] = glprim->get_opacity();

    glprim->update_glbuf(); // in case data has changed

    if (m_change_prog)
      glprim->update_shader(m_viewmode);

    glprim->get_program()->bind();
    m_ubo.bind();

    int offset = 0; // TODO: replace with one write
    offset = m_ubo.write(offset, m_udata.mvpmtx.data(), sizeof(Matrix4f));
    offset = m_ubo.write(offset, m_udata.vpmtx.data(), sizeof(Matrix4f));
    offset = m_ubo.write(offset, m_udata.modelmtx.data(), sizeof(Matrix4f));
    offset = m_ubo.write(offset, m_udata.normalmtx.data(), sizeof(Matrix4f));
    offset = m_ubo.write(offset, m_udata.vinvmtx.data(), sizeof(Matrix4f));
    offset = m_ubo.write(offset, m_udata.eyepos.data(), sizeof(Vector4f));

    offset = m_ubo.write(offset, m_udata.ambient.data(), sizeof(Vector4f));
    offset = m_ubo.write(offset, m_udata.diffuse.data(), sizeof(Vector4f));
    offset = m_ubo.write(offset, m_udata.specular.data(), sizeof(Vector4f));
    offset = m_ubo.write(offset, m_udata.options.data(), sizeof(Vector4f));
    
    Vector4f l1 = m_udata.vinvmtx * Vector4f(0.0, 0.0, 10.0, 1.0);
    glprim->get_program()->setUniformValue("lights[0].pos", QVector4D(l1[0], l1[1], l1[2], l1[3]));
    glprim->get_program()->setUniformValue("lights[0].col", QVector4D(0.2, 0.9, 0.9, 1.0));
    glprim->get_program()->setUniformValue("lights[1].pos", QVector4D(0.0, 0.0, 0.0, 0.0));
    glprim->get_program()->setUniformValue("lights[1].col", QVector4D(0.0, 0.0, 0.0, 0.0));

    glprim->get_vao().bind();

    reset_viewmode();
    if (m_viewmode == ShaderManager::PARTICLE)
    {
      glprim->get_program()->setUniformValue("pt_scale", float(14*window_dim()[1]*m_near));
      if (glprim->is_pointcloud())
      {
        glprim->get_program()->setUniformValue(
        "pt_radius",
         GLfloat(static_cast<GLPointCloud *>(glprim)->get_pointcloud()->get_radius()));
      }
      else
        glprim->get_program()->setUniformValue("pt_radius", 0.02f);
      glDrawArrays(GL_POINTS, 0, glprim->get_num_vertices());
    }
    else
    {
      glDrawElements(GL_TRIANGLES, glprim->get_num_indices(), GL_UNSIGNED_INT, 0);
    }

    glprim->get_vao().release();
    m_ubo.release();
    glprim->get_program()->release();
  }
  m_change_prog = false;
}

void SimWindow::keyPressEvent(QKeyEvent *event)
{
  int key = event->key();
  m_context->makeCurrent(this);
  switch (key)
  {
    case Qt::Key_H:
      toggle_shortcuts();
      break;
    case Qt::Key_D:
      make_dynamic();
      break;
    case Qt::Key_C:
      make_static();
      break;
    case Qt::Key_W:
      change_viewmode(ViewMode::WIREFRAME);
      break;
    case Qt::Key_S:
      change_viewmode(ViewMode::PHONG);
      break;
    case Qt::Key_P:
      change_viewmode(ViewMode::PARTICLE);
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

