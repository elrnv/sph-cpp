#ifndef SIMWINDOW_H
#define SIMWINDOW_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include <thread>
#include "openglwindow.h"
#include "uniformbuffer.h"
#include "shadermanager.h"
#include "dynamicsmanager.h"
#include "geometrymanager.h"
#include "materialmanager.h"
#include "glprimitive.h"

// forward declarations
class SPHGrid;

typedef std::vector< GLPrimitivePtr > GLPrimPtrVec;

class SimWindow : public OpenGLWindow
{
public:
  SimWindow();
  ~SimWindow();

  void init();
  void render();

  void load_model(int i);
  void init_bbox();
  void draw_wire_bbox();

  typedef ShaderManager::ShaderType ViewMode;
  void change_viewmode(ViewMode vm);
  void reset_viewmode();
  void toggle_dynamics();
  void toggle_simulation();
  void toggle_halos();
  void toggle_bbox();
  void clear_cache();

  void clear_dynamics();

  inline void toggle_shortcuts();

public slots:
  void onClose();

protected:
  void keyPressEvent(QKeyEvent *event);

private:
  bool m_show_shortcuts;
  bool m_dynamics;
  bool m_show_bbox;
  bool m_show_halos;

  UniformBuffer m_ubo;

  struct UBOData
  {
    Matrix4f mvpmtx;
    Matrix4f vpmtx;
    Matrix4f modelmtx;
    Matrix4f normalmtx; // only linear part is used
    Matrix4f vinvmtx;
    Vector4f eyepos;
  } m_udata;

  GLPrimPtrVec m_glprims;

  ViewMode m_viewmode;
  bool m_change_prog;   // flag true if m_viewmode is recently changed

  ShaderManager   m_shaderman;
  DynamicsManager m_dynman;
  GeometryManager m_geoman;
  MaterialManager m_matman;

  // bounding box
  QOpenGLVertexArrayObject m_bbox_vao;
  QOpenGLShaderProgram    *m_bbox_prog;

  std::thread m_sim_thread;
};

#endif // SIMWINDOW_H
