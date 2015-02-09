#include <QDebug>
#include <QtOpenGL/qgl.h>
#include <QtGui/QFont>
#include <QtGui/QFontDatabase>
#include <QtGui/QPainter>
#include <QtGui/QImage>
#include <QtGui/QOpenGLBuffer>
#include "eigen.h"
#include "gltext.h"

// This source must be at the head of compilation chain
namespace gl
{
  GLTextBuffer topleft;
  GLTextBuffer topright;
  GLTextBuffer bottomleft;
  GLTextBuffer bottomright;
  bool enable_text = false;
};

static uint32_t
nearest_pow2(uint32_t num)
{
  uint32_t n = num > 0 ? num - 1 : 0;

  n |= n >> 1;
  n |= n >> 2;
  n |= n >> 4;
  n |= n >> 8;
  n |= n >> 16;
  n++;

  return n;
}

// GLTextPainter stuff
GLTextPainter::GLTextPainter()
  : m_prog(this)
  , m_chars(END_CHAR - START_CHAR + 1)
{ }

void GLTextPainter::init()
{
	initializeOpenGLFunctions();

  m_prog.addShaderFromSourceFile(QOpenGLShader::Vertex,   ":/text.vert");
  m_prog.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/text.frag");
  m_prog.link();

  int id = QFontDatabase::addApplicationFont(":/proggyclean.ttf");
  QString family = QFontDatabase::applicationFontFamilies(id).at(0);
  QFont font(family, 16);
  font.setStyleHint(QFont::Monospace, QFont::NoAntialias);
  
  QFontMetrics metric(font);
  m_height = metric.height();

  int height_pow2 = nearest_pow2(m_height);

  QPainter painter;

  for (char c = START_CHAR; c <= END_CHAR; ++c)
  {
    QChar ch(c);
    FontChar fc;

    int width = metric.width(c);
    int width_pow2 = nearest_pow2(width);

    GLfloat s0 = 0.0f;
    GLfloat s1 = GLfloat(width) / GLfloat(width_pow2);

    GLfloat t0 = 0.0f;
    GLfloat t1 = -GLfloat(m_height) / GLfloat(height_pow2);

    fc.width = width;

    // paint char to QImage
    QImage final_img(width_pow2, height_pow2, QImage::Format_ARGB32);
    final_img.fill(Qt::transparent);

    painter.begin(&final_img);
    painter.setRenderHints(QPainter::Antialiasing | QPainter::HighQualityAntialiasing | QPainter::TextAntialiasing, false);
    painter.setFont(font);
    painter.setPen(Qt::black);
    painter.drawText(0, metric.ascent(), QString(c));
    painter.end();

    // create OpenGL texture
    glGenTextures(1, &fc.tex_id);
    glBindTexture(GL_TEXTURE_2D, fc.tex_id);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    final_img = QGLWidget::convertToGLFormat(final_img);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, final_img.width(), final_img.height(), 
      0, GL_RGBA, GL_UNSIGNED_BYTE, final_img.bits());

    GLfloat vtx_data[16] = {
      /* x, y, u, v */
      0,              0,                 s0, t0, // top-left
      GLfloat(width), 0,                 s1, t0, // top-right
      0,              GLfloat(m_height), s0, t1, // bottom-left
      GLfloat(width), GLfloat(m_height), s1, t1  // bottom-right
    };

    if (!m_width_vao_map.count(width))
    {
      QOpenGLVertexArrayObject *vao = new QOpenGLVertexArrayObject(this);
      m_width_vao_map.insert(std::pair<int, QOpenGLVertexArrayObject*>(width, vao));
      vao->create();
    }

    fc.vao = m_width_vao_map[width];

    fc.vao->bind();

    QOpenGLBuffer vtxbuf(QOpenGLBuffer::VertexBuffer);
    vtxbuf.create();
    vtxbuf.setUsagePattern( QOpenGLBuffer::StaticDraw );
    vtxbuf.bind();
    vtxbuf.allocate( vtx_data, sizeof( vtx_data ) );

    m_prog.enableAttributeArray( "vtx" );
    m_prog.setAttributeBuffer( "vtx", GL_FLOAT, 0, 4 );

    GLubyte idx_data[6] = { 2, 1, 0, 1, 2, 3 };
    QOpenGLBuffer idxbuf(QOpenGLBuffer::IndexBuffer);
    idxbuf.create();
    idxbuf.setUsagePattern( QOpenGLBuffer::StaticDraw );
    idxbuf.bind();
    idxbuf.allocate( idx_data, sizeof( idx_data ) );

    fc.vao->release();
    m_chars[c - START_CHAR] = fc;
  }
}

GLTextPainter::~GLTextPainter()
{
  for ( auto &x : m_width_vao_map )
    delete x.second;
}

void GLTextPainter::set_screen_size(const Vector2f &dim)
{
  m_prog.bind();
  m_prog.setUniformValue("scale", QVector2D(2.0/dim[0], -2.0/dim[1]));
  m_prog.release();
  m_screen_size = dim;
}

void GLTextPainter::draw_text()
{
  glEnable(GL_BLEND);
  glDisable(GL_DEPTH_TEST);
  glDepthMask(GL_TRUE);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  m_prog.bind();

  GLfloat left(2.0f);
  GLfloat right(m_screen_size[0] - 2.0f);
  GLfloat top(2.0f);
  GLfloat bottom(m_screen_size[1] - 2.0f - m_height);

  float opacity = 0.7f;

  // print top left buffer
  GLfloat x = left;
  GLfloat y = top;
  for ( gl::TextBlock &tb : gl::topleft )
  {
    glActiveTexture(0);
    m_prog.setUniformValue("color", QVector4D(tb.color[0], tb.color[1], tb.color[2], opacity));
    m_prog.setUniformValue("ypos", y);

    unsigned int text_len = tb.text.length();

    for (unsigned int i = 0; i < text_len; ++i)
    {
      char c = tb.text[i];

      if (c == '\n')
      {
        y += m_height;
        m_prog.setUniformValue("ypos", y);
        x = left;
        continue;
      }

      if (c < START_CHAR || c > END_CHAR)
        c = ' ';

      m_prog.setUniformValue("xpos", x);
      FontChar fc = m_chars[c - START_CHAR];
      glBindTexture(GL_TEXTURE_2D, fc.tex_id);
       
      fc.vao->bind();

      glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_BYTE, 0);

      fc.vao->release();

      x += fc.width;
    }
  }

  // print top right buffer
  x = right;
  y = top;

  for ( gl::TextBlock &tb : gl::topright )
  {
    glActiveTexture(0);
    m_prog.setUniformValue("color", QVector4D(tb.color[0], tb.color[1], tb.color[2], opacity));
    m_prog.setUniformValue("ypos", y);

    unsigned int text_len = tb.text.length();

    for (unsigned int j = text_len; j > 0; --j)
    {
      unsigned int i = j-1;
      char c = tb.text[i];

      if (c == '\n')
      {
        y += m_height;
        m_prog.setUniformValue("ypos", y);
        x = right;
        continue;
      }

      if (c < START_CHAR || c > END_CHAR)
        c = ' ';

      FontChar fc = m_chars[c - START_CHAR];
      m_prog.setUniformValue("xpos", x-fc.width);
      glBindTexture(GL_TEXTURE_2D, fc.tex_id);

      fc.vao->bind();

      glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_BYTE, 0);

      fc.vao->release();

      x -= fc.width;
    }
  }
  
  // print bottom right buffer
  x = right;
  y = bottom;

  for ( gl::TextBlock &tb : gl::bottomright )
  {
    glActiveTexture(0);
    m_prog.setUniformValue("color", QVector4D(tb.color[0], tb.color[1], tb.color[2], opacity));
    m_prog.setUniformValue("ypos", y);

    unsigned int text_len = tb.text.length();

    for (unsigned int j = text_len; j > 0; --j)
    {
      unsigned int i = j-1;
      char c = tb.text[i];

      if (c == '\n')
      {
        y -= m_height;
        m_prog.setUniformValue("ypos", y);
        x = right;
        continue;
      }

      if (c < START_CHAR || c > END_CHAR)
        c = ' ';

      FontChar fc = m_chars[c - START_CHAR];
      m_prog.setUniformValue("xpos", x-fc.width);
      glBindTexture(GL_TEXTURE_2D, fc.tex_id);

      fc.vao->bind();

      glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_BYTE, 0);

      fc.vao->release();

      x -= fc.width;
    } // for each char
  }

  // print bottom left buffer
  x = left;
  y = bottom;

  for ( gl::TextBlock &tb : gl::bottomleft )
  {
    glActiveTexture(0);
    m_prog.setUniformValue("color", QVector4D(tb.color[0], tb.color[1], tb.color[2], opacity));
    m_prog.setUniformValue("ypos", y);

    unsigned int text_len = tb.text.length();

    for (unsigned int i = 0; i < text_len; ++i)
    {
      char c = tb.text[i];

      if (c == '\n')
      {
        y -= m_height;
        m_prog.setUniformValue("ypos", y);
        x = left;
        continue;
      }

      if (c < START_CHAR || c > END_CHAR)
        c = ' ';

      m_prog.setUniformValue("xpos", x);
      FontChar fc = m_chars[c - START_CHAR];
      glBindTexture(GL_TEXTURE_2D, fc.tex_id);
       
      fc.vao->bind();

      glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_BYTE, 0);

      fc.vao->release();

      x += fc.width;
    }
  }
  m_prog.release();
  glDisable(GL_BLEND);
  glEnable(GL_DEPTH_TEST);

}
