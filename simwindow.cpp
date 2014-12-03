#include <QtGui/QScreen>
#include <QtGui/QKeyEvent>
#include <string>
#include "util.h"
#include "mesh.h"
#include "simwindow.h"
#include "glpointcloud.h"
#include "sphgrid.h"
#include "types.h"

#define STRINGIZE_VALUE(x) #x
#define STRINGIZE(x) STRINGIZE_VALUE(x)

SimWindow::SimWindow()
  : m_show_shortcuts(true) // immediately toggled below
  , m_dynamics(false)
  , m_show_bbox(true)
  , m_grid(NULL)
  , m_viewmode(ViewMode::ADDITIVE_PARTICLE)
  , m_change_prog(true)
  , m_shaderman(this)
  , m_bbox_vao(this)
{
  toggle_shortcuts();
}

SimWindow::~SimWindow()
{
}

void 
SimWindow::onClose()
{
  clear_dynamics();
}

void 
SimWindow::clear_dynamics()
{
  m_dynamics = false;
  set_animating(false);
  if (!m_grid)
    return;

  m_dynman.un_pause();

  m_dynman.request_stop(); // stop thread
  if (m_sim_thread.joinable())
    m_sim_thread.join();
  m_dynman.unrequest_stop();
  delete m_grid; m_grid = NULL;
}

void 
SimWindow::init()
{
  OpenGLWindow::init();
  m_shaderman.init();

	m_ubo.create();
	m_ubo.setUsagePattern( UniformBuffer::DynamicDraw );
	m_ubo.bind();
  m_ubo.allocate( sizeof(m_udata) );
  m_ubo.bindToIndex();

  load_model(0);

  init_bbox();
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
}

void 
SimWindow::load_model(int i)
{
  clear_dynamics();
  m_dynman.clear();
  m_geoman.clear();
  std::string cfg_filename = 
    std::string(STRINGIZE(CONFIGDIR)) + "/scene" + std::to_string(i) + ".cfg";
  bool loaded = Util::loadScene(cfg_filename, m_matman, m_geoman, m_dynman);

  if (!loaded) // nothing loaded, nothing to do
    return;

  if (global::sceneset.normalize)
  {
    m_dynman.normalize_models();
    m_dynman.transform_models( Affine3f::Identity() *
        AngleAxisf(global::sceneset.rotx, Vector3f::UnitX()));
    m_dynman.transform_models( Affine3f::Identity() *
        AngleAxisf(global::sceneset.roty, Vector3f::UnitY()));

    m_dynman.normalize_models(
        global::sceneset.padx,
        global::sceneset.pady,
        global::sceneset.padz);

    m_geoman.normalize_models();
    m_geoman.transform_models(Affine3f::Identity() *
        AngleAxisf(global::sceneset.rotx, Vector3f::UnitX()));
    m_geoman.transform_models(Affine3f::Identity() *
        AngleAxisf(global::sceneset.roty, Vector3f::UnitY()));

    m_geoman.normalize_models(
        global::sceneset.padx,
        global::sceneset.pady,
        global::sceneset.padz);
  }

  m_dynman.init_fluids(UnitBox);
  m_dynman.generate_fluiddatas();
  m_dynman.prepare_vis_data();

  glclear_tr(); // clear dynamics text buffer
  m_dynman.glprint_fluids(m_matman);

  m_udata.modelmtx.setIdentity();
  m_udata.normalmtx.block(0,0,3,3) = m_udata.modelmtx.block(0,0,3,3).inverse().transpose();

  m_glprims.clear();
  Size num_prims = 
      m_dynman.get_num_fluids() 
    + m_geoman.get_num_meshes()
    + m_geoman.get_num_pointclouds();
  m_glprims.reserve( num_prims );
  Util::loadGLData( m_glprims, m_ubo, m_shaderman, m_matman, m_geoman, m_dynman );
  change_viewmode(m_viewmode);
}

void 
SimWindow::change_viewmode(ViewMode vm)
{
  m_viewmode = vm;
  m_change_prog = true;
  reset_viewmode();
}

void 
SimWindow::reset_viewmode()
{
  glEnable(GL_BLEND);

  if (m_viewmode == ViewMode::ADDITIVE_PARTICLE) 
  {
    glEnable(GL_DEPTH_TEST);
    //glDepthFunc(GL_NEVER);
    glDepthMask(GL_FALSE);
	  glDisable(GL_CULL_FACE);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glEnable(GL_PROGRAM_POINT_SIZE);
  }
  else if (m_viewmode == ViewMode::PARTICLE) 
  {
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
	  glDisable(GL_CULL_FACE);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_PROGRAM_POINT_SIZE);
  }
  else
  {
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glEnable(GL_CULL_FACE);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_PROGRAM_POINT_SIZE);
  }
}

void 
SimWindow::toggle_bbox()
{
  m_show_bbox = !m_show_bbox;
  renderLater();
}

void 
SimWindow::toggle_halos()
{
  for ( auto glprim : m_glprims )
  {
    if (!glprim->is_pointcloud())
      continue;

    GLPointCloudPtr glpc = boost::static_pointer_cast<GLPointCloud>(glprim);
    glpc->toggle_halos();
  }
}

void 
SimWindow::clear_cache()
{
  m_dynman.clear_cache();
  glprintf_trc(CYAN, "Cache cleared\n");
}

void 
SimWindow::toggle_dynamics()
{
  bool temp = m_dynamics;

  clear_dynamics();

  m_dynamics = !temp; // toggle

  if (m_dynamics)
  {
    // Create simulation grid
    m_grid = new SPHGrid(UnitBox, m_dynman);
    m_grid->init();

    if (!m_dynman.get_bounds().size())
      m_dynman.add_default_boundary(*m_grid); // always need a boundary

    // run simulation
    m_sim_thread = std::thread(&DynamicsManager::run, &m_dynman, m_grid);
  }
  set_animating(m_dynamics);
}

void 
SimWindow::toggle_simulation()
{ 
  m_dynman.toggle_pause();
}

void 
SimWindow::render()
{
  glDepthMask(GL_TRUE); // depthmask needs to be true before clearing the depth bit
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

  Affine3f vtrans(get_view_trans()); // view transformation

  AffineCompact3f &mvtrans = get_view_trans();
  if (m_viewmode == ViewMode::PARTICLE)
  { // explicitly transform and sort by depth
    for ( auto glprim : m_glprims )
      glprim->update_glbuf_withsort(mvtrans, mvtrans);
  //  vtrans.setIdentity(); // reset the view matrix
  }
  else
  {
    for ( auto glprim : m_glprims )
      glprim->update_glbuf_nosort();
  }

  // sort primitives (why?)
  /*
  if (m_viewmode != ViewMode::ADDITIVE_PARTICLE)
  {
    std::sort(m_glprims.begin(), m_glprims.end(),
        [mvtrans](const GLPrimitive *p1, const GLPrimitive *p2)
        { return (mvtrans * p1->get_closest_pt())[2] < (mvtrans * p2->get_closest_pt())[2]; });
  }*/

  // write uniform buffer
  m_udata.vpmtx = (get_proj_trans() * vtrans).matrix();
  m_udata.mvpmtx = m_udata.vpmtx; // ignoring unused model transformations
  m_udata.vinvmtx = vtrans.inverse().matrix();
  m_udata.eyepos = m_udata.vinvmtx * Vector4f(0.0f, 0.0f, 0.0f, 1.0f);

  m_ubo.bind();
  m_ubo.write(0, &m_udata, sizeof( m_udata )); // write uniform buffer object
  m_ubo.release();

  glFinish(); // Finish writing uniform buffer before drawing

  if (m_show_bbox)
    draw_wire_bbox();

  reset_viewmode();

  for ( auto glprim : m_glprims )
  {
    if (m_change_prog)
      glprim->update_shader(m_viewmode, m_shaderman);

    glprim->get_program()->bind();

    glprim->get_program()->setUniformValue("ambient",  glprim->get_ambient());
    glprim->get_program()->setUniformValue("diffuse",  glprim->get_diffuse());
    glprim->get_program()->setUniformValue("specular", glprim->get_specular());
    glprim->get_program()->setUniformValue("specpow",  glprim->get_specpow());
    glprim->get_program()->setUniformValue("opacity",  glprim->get_opacity());

    Vector4f l1 = m_udata.vinvmtx * Vector4f(0.0, 0.0, 10.0, 1.0);
    glprim->get_program()->setUniformValue("lights[0].pos", QVector4D(l1[0], l1[1], l1[2], l1[3]));
    glprim->get_program()->setUniformValue("lights[0].col", QVector4D(0.9, 0.9, 0.9, 1.0));
    glprim->get_program()->setUniformValue("lights[1].pos", QVector4D(0.0, 0.0, 0.0, 0.0));
    glprim->get_program()->setUniformValue("lights[1].col", QVector4D(0.0, 0.0, 0.0, 0.0));

    glprim->get_vao().bind();
    if (m_viewmode == ViewMode::PARTICLE ||
        m_viewmode == ViewMode::ADDITIVE_PARTICLE)
    {
      glprim->get_program()->setUniformValue("pt_scale", float(14.5*window_dim()[1]*m_near));
      if (glprim->is_pointcloud())
      {
        GLPointCloudPtr glpc = boost::static_pointer_cast<GLPointCloud>(glprim);
        glprim->get_program()->setUniformValue( "pt_radius", GLfloat(glpc->get_radius()));
        if (!glpc->is_halos())
          glprim->get_program()->setUniformValue( "pt_halo", GLfloat(glpc->get_radius()));
        else
          glprim->get_program()->setUniformValue("pt_halo", GLfloat(glpc->get_halo_radius()));
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


void 
SimWindow::toggle_shortcuts()
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
  glprintf_blc(BLUE,  "    Q - particles\n");
  glprintf_blc(BLUE,  "    A - additive particles\n");
  glprintf_bl("  Scenes: \n");
  glprintf_blc(BLUE,  "  0-9 - Custom scene configurations\n");
  glprintf_bl("  Dynamics: \n");
  glprintf_blc(RED,  "    D - start/stop\n");
  glprintf_blc(RED,  "    G - pause/resume\n");
  glprintf_blc(RED,  "    C - clear cache\n");
  glprintf_blc(CYAN, "    H - show/hide halos\n");
  glprintf_blc(CYAN, "    B - show/hide bounding box\n");
}

void 
SimWindow::keyPressEvent(QKeyEvent *event)
{
  int key = event->key();
  //m_context->makeCurrent(this);
  switch (key)
  {
    case Qt::Key_T:
      toggle_shortcuts();
      break;
    case Qt::Key_D:
      toggle_dynamics();
      break;
    case Qt::Key_G:
      toggle_simulation();
      break;
    case Qt::Key_H:
      toggle_halos();
      break;
    case Qt::Key_C:
      clear_cache();
      break;
    case Qt::Key_W:
      change_viewmode(ViewMode::NORMALS);
      break;
    case Qt::Key_S:
      change_viewmode(ViewMode::PHONG);
      break;
    case Qt::Key_Q:
      change_viewmode(ViewMode::PARTICLE);
      break;
    case Qt::Key_A:
      change_viewmode(ViewMode::ADDITIVE_PARTICLE);
      break;
    case Qt::Key_B:
      toggle_bbox();
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

void 
SimWindow::init_bbox()
{
  GLfloat vertices[24] = {
    -1, -1, -1,   1, -1, -1,
    -1, -1,  1,   1, -1,  1,
    -1,  1,  1,   1,  1,  1,
    -1,  1, -1,   1,  1, -1
  };

  GLubyte indices[24] = {
    0, 2, 2, 4, 4, 6, 6, 0,
    1, 3, 3, 5, 5, 7, 7, 1,
    0, 1, 2, 3, 4, 5, 6, 7
  };

  m_bbox_vao.create();
  m_bbox_vao.bind();

  QOpenGLBuffer vtxbuf(QOpenGLBuffer::VertexBuffer);
  vtxbuf.create();
  vtxbuf.setUsagePattern( QOpenGLBuffer::StaticDraw );
  vtxbuf.bind();
  vtxbuf.allocate( vertices, sizeof( vertices ) );

  m_bbox_prog = m_shaderman.get_flat_shader();
  m_bbox_prog->enableAttributeArray( "pos" );
  m_bbox_prog->setAttributeBuffer( "pos", GL_FLOAT, 0, 3 );

  QOpenGLBuffer idxbuf(QOpenGLBuffer::IndexBuffer);
  idxbuf.create();
  idxbuf.setUsagePattern( QOpenGLBuffer::StaticDraw );
  idxbuf.bind();
  idxbuf.allocate( indices, sizeof( indices ) );

  m_bbox_vao.release();

  m_ubo.bindToProg(m_bbox_prog->programId(), "Globals");
}

void 
SimWindow::draw_wire_bbox()
{
  glEnable(GL_BLEND);
  glEnable(GL_DEPTH_TEST);
  glDepthMask(GL_TRUE);
  glDisable(GL_CULL_FACE);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDisable(GL_PROGRAM_POINT_SIZE);

  m_bbox_prog->bind();
  m_bbox_prog->setUniformValue("diffuse", QVector3D(0.3f, 0.4f, 0.5f));
  m_bbox_prog->setUniformValue("opacity", 0.5f);

  m_bbox_vao.bind();
  glDrawElements(GL_LINES, 24, GL_UNSIGNED_BYTE, 0);
  m_bbox_vao.release();

  m_bbox_prog->release();
}

