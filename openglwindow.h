#ifndef OPENGLWINDOW_H
#define OPENGLWINDOW_H

#include <QDebug>
#include <QtGui/QWindow>
#include <QtGui/QOpenGLFunctions_3_3_Core>
#include <QtGui/QMatrix4x4>
#include <QTime>
#include "eigen.h"
#include "gltext.h"

class OpenGLWindow : public QWindow, protected QOpenGLFunctions_3_3_Core
{
  Q_OBJECT

public:
  explicit OpenGLWindow(QWindow *parent = 0);
  ~OpenGLWindow();

  virtual void init();
  virtual void reshape();
  virtual void render();
  void set_animating(bool animating);

public slots:
  void renderLater();
  void renderNow();

protected:
  Vector2f window_dim();

  bool event(QEvent *event);

  void resizeEvent(QResizeEvent *event);
  void exposeEvent(QExposeEvent *event);
  void showEvent(QShowEvent *event);

  void mousePressEvent(QMouseEvent *event);
  void mouseReleaseEvent(QMouseEvent *event);
  void mouseMoveEvent(QMouseEvent *event);
  void wheelEvent(QWheelEvent *event);
  void keyPressEvent(QKeyEvent *event);

  Projective3f &get_proj_trans() { return m_P; }
  AffineCompact3f &get_view_trans() { return m_V; }

  void reset_view();
  void recompute_view();

  float m_near; // near plane
  float m_far;  // far plane
  float m_fov;  // field of view

  unsigned int m_frame;
  QTime        m_time;

private:

  bool m_update_pending;
  bool m_animating;
  bool m_rotation_control;

  int m_prev_x, m_prev_y; // previous cursor position at moust click

  float m_hra; // horizontal rotation angle
  float m_vra; // vertical rotation angle
  float m_zoom; // zoom level

  Projective3f m_P; // perspective projection
  AffineCompact3f m_V; // view transform

  QOpenGLContext *m_context;
};

#endif // OPENGLWINDOW_H
