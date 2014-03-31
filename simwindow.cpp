#include <QtGui/QScreen>
#include <QtGui/QKeyEvent>
#include <string>
#include "scene.h"
#include "util.h"
#include "mesh.h"
#include "simwindow.h"
#include "fluid.h"

void SimWindow::toggle_shortcuts()
{
  m_show_shortcuts = !m_show_shortcuts;
  glclear_bl();

  if (!m_show_shortcuts)
  {
    glprintf_blc(BLUE, "T - show/hide shortcuts\n");
    return;
  }
  glprintf_bl("Shortcuts:\n");
  glprintf_blc(BLUE, "  Y - enable/disable text\n");
  glprintf_blc(BLUE, "  T - show/hide shortcuts\n");
  glprintf_blc(BLUE, "  R - reset view\n");
  glprintf_bl("  View Modes: \n");
  glprintf_blc(GREEN, "    W - wireframe\n");
  glprintf_blc(RED,   "    S - phong\n");
  glprintf_blc(BLUE,  "    P - particles\n");
  glprintf_bl("  Models: \n");
  glprintf_blc(GREEN, "    1 - cube\n");
  glprintf_blc(GREEN, "    2 - cow\n");
  glprintf_blc(GREEN, "    3 - bunny\n");
  glprintf_blc(BLUE,  "    4 - sparse point cloud\n");
  glprintf_blc(BLUE,  "    5 - normal point cloud\n");
  glprintf_blc(BLUE,  "    6 - dense point cloud\n");
  glprintf_bl("  Dynamics: \n");
  glprintf_blc(RED,  "    D - start\n");
  glprintf_blc(RED,  "    C - stop\n");
  glprintf_blc(CYAN, "    H - show/hide halos\n");
}

SimWindow::SimWindow()
  : m_show_shortcuts(true) // immediately toggled below
  , m_grid(NULL)
  , m_viewmode(ShaderManager::PARTICLE)
  , m_change_prog(true)
  , m_shaderman(this)
{
  toggle_shortcuts();
}

SimWindow::~SimWindow()
{
  clear_dynamics();
}

void SimWindow::clear_dynamics()
{
  if (!m_grid)
    return;

  m_grid->request_stop(); // stop thread
  m_sim_thread.join();
  delete m_grid; m_grid = 0;
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

  load_model(0);
}

void SimWindow::load_model(int i)
{
  clear_dynamics();
  SceneNode *scene = Util::loadScene("data/scene" + std::to_string(i) + ".cfg");

  if (!scene)
    return;

  scene->normalize_model();
  scene->rotate(global::sceneset.rotx, Vector3f::UnitX());
  scene->rotate(global::sceneset.roty, Vector3f::UnitY());
  scene->flatten();

  scene->normalize_model(
      global::sceneset.padx,
      global::sceneset.pady,
      global::sceneset.padz);
  scene->flatten();

  scene->cube_bbox();
  m_udata.modelmtx.setIdentity();

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
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

  if (m_viewmode == ViewMode::PARTICLE) 
  {
    glDisable(GL_DEPTH_TEST);
    //glDepthFunc(GL_NEVER);
	  glDisable(GL_CULL_FACE);
    glBlendFunc(GL_ONE, GL_ONE);
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

void SimWindow::stop_dynamics()
{
  clear_dynamics();
  set_animating(false);
}

void SimWindow::toggle_halos()
{
  for ( const GLPrimitivePtr &prim_ptr : m_glprims )
  {
    GLPrimitive *glprim = prim_ptr.get();
    if (!glprim->is_pointcloud())
      continue;

    GLPointCloud *glpc = static_cast<GLPointCloud *>(glprim);
    glpc->toggle_halos();
  }
}

void SimWindow::start_dynamics()
{
  clear_dynamics();
  glclear_tr(); // clear dynamics text buffer

  // Create simulation grid
  m_grid = new UniformGrid(Vector3f(-1,-1,-1), Vector3f(1,1,1));

  for ( const GLPrimitivePtr &prim_ptr : m_glprims )
  {
    GLPrimitive *glprim = prim_ptr.get();
    if (!glprim->is_pointcloud())
      continue;

    GLPointCloud *glpc = static_cast<GLPointCloud *>(glprim);
    if (glpc->is_dynamic())
      m_grid->add_fluid(glpc);
  }

  // run simulation
  m_sim_thread = std::thread(&UniformGrid::run, m_grid);

  set_animating(true);
}

void SimWindow::render()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

  reset_viewmode();
  int i = 0;

  Affine3f vtrans(get_view_trans()); // view transformation
  m_udata.vpmtx = (get_proj_trans() * vtrans).matrix();
  m_udata.mvpmtx = m_udata.vpmtx * m_udata.modelmtx;
  m_udata.vinvmtx = vtrans.inverse().matrix();
  m_udata.eyepos = m_udata.vinvmtx * Vector4f(0.0f, 0.0f, 0.0f, 1.0f);

  m_udata.normalmtx.block(0,0,3,3) = m_udata.modelmtx.block(0,0,3,3).inverse().transpose();

  m_ubo.bind();
  m_ubo.write(0, &m_udata, sizeof( m_udata )); // write uniform buffer object
  m_ubo.release();

  glFinish(); // Finish writing uniform buffer before drawing

#if 0 // sorting
  AffineCompact3f mvtrans = AffineCompact3f(vtrans.matrix() * m_udata.modelmtx);

  for ( const GLPrimitivePtr &prim_ptr : m_glprims )
  {
    GLPrimitive *glprim = prim_ptr.get();

    if (m_viewmode == ShaderManager::PARTICLE)
      glprim->sort_by_depth(mvtrans);

    ++i;
  }

  std::sort(m_glprims.begin(), m_glprims.end(),
      [mvtrans](const GLPrimitivePtr &p1, const GLPrimitivePtr &p2)
      { return (mvtrans * p1->get_closest_pt())[2] < (mvtrans * p2->get_closest_pt())[2]; });
#endif

  for ( const GLPrimitivePtr &prim_ptr : m_glprims )
  {
    GLPrimitive *glprim = prim_ptr.get();

    glprim->update_glbuf(); // in case data has changed

    if (m_change_prog)
      glprim->update_shader(m_viewmode);

    glprim->get_program()->bind();

    glprim->get_program()->setUniformValue("ambient", glprim->get_ambient());
    glprim->get_program()->setUniformValue("diffuse", glprim->get_diffuse());
    glprim->get_program()->setUniformValue("specular", glprim->get_specular());
    glprim->get_program()->setUniformValue("specpow", glprim->get_specpow());
    glprim->get_program()->setUniformValue("opacity", glprim->get_opacity());

    Vector4f l1 = m_udata.vinvmtx * Vector4f(0.0, 0.0, 10.0, 1.0);
    glprim->get_program()->setUniformValue("lights[0].pos", QVector4D(l1[0], l1[1], l1[2], l1[3]));
    glprim->get_program()->setUniformValue("lights[0].col", QVector4D(0.9, 0.9, 0.9, 1.0));
    glprim->get_program()->setUniformValue("lights[1].pos", QVector4D(0.0, 0.0, 0.0, 0.0));
    glprim->get_program()->setUniformValue("lights[1].col", QVector4D(0.0, 0.0, 0.0, 0.0));

    glprim->get_vao().bind();
    if (m_viewmode == ShaderManager::PARTICLE)
    {
      glprim->get_program()->setUniformValue("pt_scale", float(14.5*window_dim()[1]*m_near));
      if (glprim->is_pointcloud())
      {
        GLPointCloud * glpc = static_cast<GLPointCloud *>(glprim);
        PointCloud *pc = glpc->get_pointcloud();
        glprim->get_program()->setUniformValue( "pt_radius", GLfloat(pc->get_radius()));
        if (pc->is_dynamic() && !glpc->is_halos())
        {
          glprim->get_program()->setUniformValue( "pt_halo", GLfloat(pc->get_radius()));
        }
        else
          glprim->get_program()->setUniformValue("pt_halo", GLfloat(pc->get_halo_radius()));
      }
      else
      {
        glprim->get_program()->setUniformValue("pt_radius", 0.02f);
        glprim->get_program()->setUniformValue("pt_halo", 0.02f);
      }

      glDrawArrays(GL_POINTS, 0, glprim->get_num_vertices());
    }
    else
    {
      glDrawElements(GL_TRIANGLES, glprim->get_num_indices(), GL_UNSIGNED_INT, 0);
    }

    glprim->get_vao().release();
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
    case Qt::Key_T:
      toggle_shortcuts();
      break;
    case Qt::Key_D:
      start_dynamics();
      break;
    case Qt::Key_H:
      toggle_halos();
      break;
    case Qt::Key_C:
      stop_dynamics();
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

