#ifndef SIMWINDOW_H
#define SIMWINDOW_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include "mesh.h"
#include "openglwindow.h"
#include "uniformbuffer.h"
#include "shadermanager.h"

typedef boost::shared_ptr<GLPrimitive> GLPrimitivePtr;
typedef std::vector< GLPrimitivePtr > GLPrimVec;

class SimWindow : public OpenGLWindow
{
public:
  SimWindow();

  void init();
  void render();

  void load_model(int i);

  typedef ShaderManager::ShaderType ViewMode;
  void update_viewmode(ViewMode vm);
  void make_dynamic();

protected:
  void keyPressEvent(QKeyEvent *event);

private:
  UniformBuffer m_global_uniform;

  struct UBOGlobals
  {
    Matrix4f mvpmtx;
    Matrix4f vpmtx;
    Matrix4f modelmtx;
    Matrix4f normalmtx;
    Matrix4f vinvmtx;
    Vector4f eyepos;
  } m_ubo;

  GLPrimVec m_glprims;

  ViewMode m_viewmode;

  ShaderManager m_shaderman;
};

#endif // SIMWINDOW_H
