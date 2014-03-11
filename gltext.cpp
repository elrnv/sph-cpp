#include <QDebug>
#include <QtOpenGL/qgl.h>
#include <QtGui/QFont>
#include <QtGui/QFontDatabase>
#include <QtGui/QPainter>
#include <QtGui/QImage>
#include <QtGui/QOpenGLBuffer>
#include "eigen.h"
#include "gltext.h"

namespace gl
{
  GLText text;
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

// GLText stuff
GLText::GLText()
  : m_prog(this)
  , m_chars(END_CHAR - START_CHAR + 1)
{ }

void GLText::init()
{
	initializeOpenGLFunctions();

  m_prog.addShaderFromSourceFile(QOpenGLShader::Vertex,   ":/text.vert");
  m_prog.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/text.frag");
  m_prog.link();

  int id = QFontDatabase::addApplicationFont(":/proggyclean.ttf");
  QString family = QFontDatabase::applicationFontFamilies(id).at(0);
  QFont font(family, 16);
  
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
      0,              0,                    s0, t0, // top-left
      GLfloat(width), 0,                    s1, t0, // top-right
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

    GLubyte idx_data[6] = { 0, 1, 2, 1, 2, 3 };
    QOpenGLBuffer idxbuf(QOpenGLBuffer::IndexBuffer);
    idxbuf.create();
    idxbuf.setUsagePattern( QOpenGLBuffer::StaticDraw );
    idxbuf.bind();
    idxbuf.allocate( idx_data, sizeof( idx_data ) );

    fc.vao->release();
    m_chars[c - START_CHAR] = fc;
  }
}

GLText::~GLText()
{
  for ( auto &x : m_width_vao_map )
    delete x.second;
}

void GLText::set_screen_size(const Vector2f &dim)
{
  m_prog.bind();
  m_prog.setUniformValue("scale", QVector2D(2.0/dim[0], -2.0/dim[1]));
  m_prog.release();
  m_screen_size = dim;
}

// printing routines
void GLText::print(const std::string &text)
{
  print(text, gl::TOP_LEFT, gl::WHITE );
}
void GLText::print(const std::string &text, gl::Corner corner)
{
  print(text, corner, gl::WHITE);
}
void GLText::print(const std::string &text, gl::Color color)
{
  print(text, gl::TOP_LEFT, color);
}
void GLText::print(const std::string &text, gl::Corner corner, gl::Color color)
{
  Vector3f c;
  switch (color)
  {
    case gl::RED:   c = Vector3f(1.0f, 0.0f, 0.0f); break;
    case gl::GREEN: c = Vector3f(0.0f, 1.0f, 0.0f); break;
    case gl::BLUE:  c = Vector3f(0.0f, 0.0f, 1.0f); break;
    case gl::BLACK: c = Vector3f(0.0f, 0.0f, 0.0f); break;
    case gl::WHITE: // FALL THROUGH
    default: c = Vector3f(1.0f, 1.0f, 1.0f); break;
  }
  switch (corner)
  {
    case gl::TOP_RIGHT: m_topright.push_back(TextBlock( text, c )); break;
    case gl::TOP_LEFT:  // FALL THROUGH
    default: m_topleft.push_back(TextBlock( text, c )); break;
  }
}

void GLText::draw_text()
{
  glEnable(GL_BLEND);
  glDisable(GL_DEPTH_TEST);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  m_prog.bind();

  GLfloat left(2.0f);
  GLfloat right(m_screen_size[0] - 2.0f);
  GLfloat top(2.0f);

  GLfloat x = left;
  GLfloat y = top;
  for ( TextBlock &tb : m_topleft )
  {
    glActiveTexture(0);
    m_prog.setUniformValue("color", QVector4D(tb.color[0], tb.color[1], tb.color[2], 1.0));
    m_prog.setUniformValue("ypos", y);

    unsigned int text_len = tb.text.length();

    for (unsigned int i = 0; i < text_len; ++i)
    {
      char c = tb.text[i];

      if (c == '\n')
      {
        y += m_height;
        m_prog.setUniformValue("ypos", y);
        x = GLfloat(2.0f);
        continue;
      }

      if (c < START_CHAR && c > END_CHAR)
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

  x = right;
  y = top;

  for ( TextBlock &tb : m_topright )
  {
    glActiveTexture(0);
    m_prog.setUniformValue("color", QVector4D(tb.color[0], tb.color[1], tb.color[2], 1.0));
    m_prog.setUniformValue("ypos", y);

    unsigned int text_len = tb.text.length();

    for (unsigned int j = text_len; j > 0; --j)
    {
      unsigned int i = j;
      char c = tb.text[i];

      if (c == '\n')
      {
        y += m_height;
        m_prog.setUniformValue("ypos", y);
        x = right;
        continue;
      }

      if (c < START_CHAR && c > END_CHAR)
        c = ' ';

      m_prog.setUniformValue("xpos", x);
      FontChar fc = m_chars[c - START_CHAR];
      glBindTexture(GL_TEXTURE_2D, fc.tex_id);

      fc.vao->bind();

      glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_BYTE, 0);

      fc.vao->release();

      x -= fc.width;
    }
  }
  m_prog.release();
  glDisable(GL_BLEND);
  glEnable(GL_DEPTH_TEST);

}
