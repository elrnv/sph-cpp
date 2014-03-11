#ifndef GLTEXT_H
#define GLTEXT_H

#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLShaderProgram>
#include <QtGui/QOpenGLFunctions_3_3_Core>
#include <boost/unordered_map.hpp>
#include "eigen.h"

namespace gl
{
  enum Corner
  {
    TOP_LEFT,
    TOP_RIGHT
  };
  enum Color
  {
    RED,
    GREEN,
    BLUE,
    BLACK,
    WHITE
  };
};

struct FontChar
{
  int width;     // width of char
  GLuint tex_id; // texture id
  QOpenGLVertexArrayObject *vao;
};

// class to manage drawing text in modern OpenGL
class GLText : public QObject, protected QOpenGLFunctions_3_3_Core
{
public:
  explicit GLText();
  ~GLText();

  void init();
  void set_screen_size(const Vector2f &dim);
  void print(const std::string &text);
  void print(const std::string &text, gl::Corner corner);
  void print(const std::string &text, gl::Color color);
  void print(const std::string &text, gl::Corner corner, gl::Color color);
  void draw_text();

protected:
  boost::unordered_map<int, QOpenGLVertexArrayObject *> m_width_vao_map; // dictionary of vaos

  QOpenGLShaderProgram m_prog;  // glsl programs used for rendering

  const static char START_CHAR = 32;
  const static char END_CHAR   = 126;
  std::vector<FontChar> m_chars; // array of characters

  float    m_height;

  struct TextBlock
  {
    TextBlock( const std::string &str, const Vector3f &c )
      : text(str), color(c){ }
    std::string text;
    Vector3f color;
  };

  std::vector<TextBlock> m_topleft;
  std::vector<TextBlock> m_topright;

  Vector2f m_screen_size;
};

namespace gl
{
  extern GLText text;
};

#endif // GLTEXT_H
