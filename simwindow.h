#ifndef SIMWINDOW_H
#define SIMWINDOW_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include "openglwindow.h"
#include "uniformbuffer.h"
#include "shadermanager.h"
#include "glprimitive.h"

typedef boost::shared_ptr<GLPrimitive> GLPrimitivePtr;
typedef std::vector< GLPrimitivePtr > GLPrimVec;

class SimWindow : public OpenGLWindow
{
public:
  SimWindow();
  ~SimWindow();

  void init();
  void render();

  void load_model(int i);

  typedef ShaderManager::ShaderType ViewMode;
  void change_viewmode(ViewMode vm);
  void reset_viewmode();
  void make_static();
  void make_dynamic();
  void toggle_halos();

  void clear_dynamics();

  inline void toggle_shortcuts();

protected:
  void keyPressEvent(QKeyEvent *event);

private:
  bool m_show_shortcuts;

  UniformBuffer m_ubo;

  struct UBOData
  {
    Matrix4f mvpmtx;
    Matrix4f vpmtx;
    Matrix4f modelmtx;
    Matrix4f normalmtx; // only linear part is used
    Matrix4f vinvmtx;
    Vector4f eyepos;

    Vector4f ambient;
    Vector4f diffuse;
    Vector4f specular;
    Vector4f options; // only the first value is used
  } m_udata;

  UnfiromGrid *m_grid; // simulation grid

  GLPrimVec m_glprims;

  ViewMode m_viewmode;
  bool m_change_prog; // flag true if m_viewmode is recently changed

  ShaderManager m_shaderman;

  std::thread m_sim_thread;
};

#endif // SIMWINDOW_H
