#include <QDebug>
#include <QtCore/QCoreApplication>
#include <QtGui/QOpenGLContext>
#include <QtGui/QMouseEvent>
#include "eigen.h"
#include "openglwindow.h"

OpenGLWindow::OpenGLWindow(QWindow *parent)
	: QWindow(parent)
  , m_near(0.1), m_far(100.0), m_fov(40.0)
  , m_frame(0)
	, m_update_pending(false)
	, m_animating(false)
	, m_rotation_control(false)
	, m_prev_x(0), m_prev_y(0)
  , m_zoom(0)
	, m_context(0)
{
	setSurfaceType(QWindow::OpenGLSurface);
  m_time.start();
}

OpenGLWindow::~OpenGLWindow()
{

}

void OpenGLWindow::init() { }


Vector2f OpenGLWindow::window_dim()
{
	const qreal retinaScale = devicePixelRatio();
	return Vector2f(width(), height()) * retinaScale;
}

void OpenGLWindow::reshape()
{
  Vector2f dim = window_dim();
	glViewport(0, 0, dim[0], dim[1]);
	m_P.setIdentity();
	m_P.perspective(m_fov, dim[0]/dim[1], m_near, m_far);
}

void OpenGLWindow::render() { }

bool OpenGLWindow::event(QEvent *event)
{
	switch (event->type()) {
		case QEvent::UpdateRequest:
			m_update_pending = false;
			renderNow();
			return true;
		default:
			return QWindow::event(event);
	}
}

void OpenGLWindow::resizeEvent(QResizeEvent *event)
{
	Q_UNUSED(event);
	reshape();
	renderNow();
}

void OpenGLWindow::exposeEvent(QExposeEvent *event)
{
	Q_UNUSED(event);
	renderNow();
}

void OpenGLWindow::showEvent(QShowEvent *event)
{
	Q_UNUSED(event);
	
	if (m_context)
		return;

	// Initialize OpenGL Context
	m_context = new QOpenGLContext(this);

	// Create platform-specific surface format
	QSurfaceFormat format;
	format.setDepthBufferSize(24);
	format.setStencilBufferSize(8);
	format.setVersion( 3, 3 );
	format.setSamples( 8 );
  format.setSwapBehavior( QSurfaceFormat::DoubleBuffer );
	format.setProfile( QSurfaceFormat::CoreProfile );

	m_context->setFormat(format);
	m_context->create();
	m_context->makeCurrent(this);

	initializeOpenGLFunctions();
	init();
  reset_view();
}

void OpenGLWindow::mousePressEvent(QMouseEvent *event)
{
	Qt::MouseButton button = event->button();
	switch (button)
	{
		case Qt::LeftButton:
			m_prev_x = event->x();
			m_prev_y = event->y();
			m_rotation_control = true;
			break;
		default:
			break;
	}
}

void OpenGLWindow::mouseReleaseEvent(QMouseEvent *event)
{
	Qt::MouseButton button = event->button();
	switch (button)
	{
		case Qt::LeftButton:
			m_rotation_control = false;
			break;
		default:
			break;
	}
}

void OpenGLWindow::mouseMoveEvent(QMouseEvent *event)
{
  int y = event->y();
  if (event->modifiers() & Qt::ShiftModifier)
  {
    m_zoom += event->y() - m_prev_y;
    m_zoom = std::max(m_zoom, -99.0f);
    m_zoom = std::min(m_zoom, 1000.0f);
    recompute_view();
    renderLater();
  }
  else if (m_rotation_control)
  {
    int x = event->x();
    m_hra = (m_hra - m_prev_x + x);
    m_hra = m_hra > 180 ? m_hra - 360 : m_hra;
    m_hra = m_hra < -180 ? m_hra + 360 : m_hra;
    m_prev_x = x;

    float vra_temp = std::max(m_vra - m_prev_y + y, -180.0f);
    m_vra = std::min(vra_temp, 180.0f);
    
    recompute_view();
    renderLater();
  }
  m_prev_y = y; // always update previous value of y (for zoom)
}

void OpenGLWindow::wheelEvent(QWheelEvent *event)
{
  QPoint num_degrees = event->angleDelta() / 8;

  if (!num_degrees.isNull())
  {
    m_zoom += 0.1*num_degrees.y();
  }

  m_zoom = std::max(m_zoom, -99.0f);
  m_zoom = std::min(m_zoom, 1000.0f);

  recompute_view();
	renderLater();
}


void OpenGLWindow::keyPressEvent(QKeyEvent *event)
{
	int key = event->key();
	switch (key)
	{
		case Qt::Key_R:
      reset_view();
			break;
		default:
			break;
	}
	
	renderLater();
}

// Immediate render trigger
void OpenGLWindow::renderNow()
{
	if (!isExposed())
		return;

	m_context->makeCurrent(this);

	render();

	m_context->swapBuffers(this);

	if (m_animating)
  {
    ++m_frame;
    if (m_frame == 100)
    {
      //fprintf(stderr, "\r%d       ", 100000 / m_time.elapsed() );
      m_time.restart();
      m_frame = 0;
    }
		renderLater();
  }
}


// Queued render trigger
void OpenGLWindow::renderLater()
{
	if (!m_update_pending) {
		m_update_pending = true;
		QCoreApplication::postEvent(this, new QEvent(QEvent::UpdateRequest));
	}
}

void OpenGLWindow::set_animating(bool animating)
{
	m_animating = animating;

	if (animating)
  {
    m_time.restart();
    m_frame = 0;
		renderLater();
  }
}

// view related functions
void OpenGLWindow::reset_view()
{ 
  m_zoom = 0;
  m_hra = 0.0f;//25.0f; // degrees
  m_vra = 0.0f;//25.0f; // degrees
  recompute_view();
}

void OpenGLWindow::recompute_view()
{
  m_V.setIdentity();
  m_V.translate(Vector3f(0,0,-5));
  m_V.scale(1 + 0.01*m_zoom);
  m_V.rotate(AngleAxisf(m_vra*RADIAN, Vector3f::UnitX()));
  m_V.rotate(AngleAxisf(m_hra*RADIAN, Vector3f::UnitY()));
}
