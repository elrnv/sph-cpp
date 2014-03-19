#ifndef GLTEXT_H
#define GLTEXT_H

#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLShaderProgram>
#include <QtGui/QOpenGLFunctions_3_3_Core>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include "eigen.h"

// class to manage drawing text in modern OpenGL
class GLTextPainter : public QObject, protected QOpenGLFunctions_3_3_Core
{
public:
  explicit GLTextPainter();
  ~GLTextPainter();

  void init();
  void set_screen_size(const Vector2f &dim);
  void draw_text(int x, int y);

protected:
  boost::unordered_map<int, QOpenGLVertexArrayObject *> m_width_vao_map; // dictionary of vaos

  QOpenGLShaderProgram m_prog;  // glsl programs used for rendering

  const static char START_CHAR = 32;
  const static char END_CHAR   = 126;

  // structure used to draw characters on screen
  struct FontChar
  {
    int width;     // width of char
    GLuint tex_id; // texture id
    QOpenGLVertexArrayObject *vao;
  };

  std::vector<FontChar> m_chars; // array of characters

  float    m_height;
  Vector2f m_screen_size;
};

namespace gl
{
  // structure used to keep track of text blocks in the buffer
  struct TextBlock
  {
    TextBlock( const std::string &str, const Vector3f &c )
      : text(str), color(c){ }
    std::string text;
    Vector3f color;
  };

  typedef std::deque<TextBlock> GLTextBuffer;

  enum Corner
  {
    TOP_LEFT,
    TOP_RIGHT,
    BOTTOM_LEFT,
    BOTTOM_RIGHT,
  };
  enum Color
  {
    RED,
    YELLOW,
    GREEN,
    CYAN,
    BLUE,
    BLACK,
    WHITE
  };

  extern GLTextBuffer topleft;
  extern GLTextBuffer topright;
  extern GLTextBuffer bottomleft;
  extern GLTextBuffer bottomright;
};

inline void _glclear(gl::Corner corner)
{
  switch (corner)
  {
    case gl::BOTTOM_RIGHT: gl::bottomright.clear(); break;
    case gl::BOTTOM_LEFT:  gl::bottomleft.clear(); break;
    case gl::TOP_RIGHT: gl::topright.clear(); break;
    case gl::TOP_LEFT:  // FALL THROUGH
    default: gl::topleft.clear(); break;
  }
}

#define glclear_tr() _glclear(gl::TOP_RIGHT)
#define glclear_tl() _glclear(gl::TOP_LEFT)
#define glclear_br() _glclear(gl::BOTTOM_RIGHT)
#define glclear_bl() _glclear(gl::BOTTOM_LEFT)

inline void _glprintf(gl::Corner corner, gl::Color color, const char *fmt, ...)
{
  char buffer[256];
  va_list args;
  va_start (args, fmt);
  vsprintf (buffer, fmt, args);
  va_end (args);

  Vector3f c;
  switch (color)
  {
    case gl::RED:   c = Vector3f(1.0f, 0.0f, 0.0f); break;
    case gl::YELLOW: c = Vector3f(1.0f, 1.0f, 0.0f); break;
    case gl::GREEN: c = Vector3f(0.0f, 1.0f, 0.0f); break;
    case gl::CYAN:  c = Vector3f(0.0f, 1.0f, 1.0f); break;
    case gl::BLUE:  c = Vector3f(0.3f, 0.4f, 1.0f); break;
    case gl::BLACK: c = Vector3f(0.0f, 0.0f, 0.0f); break;
    case gl::WHITE: // FALL THROUGH
    default: c = Vector3f(1.0f, 1.0f, 1.0f); break;
  }
  std::deque<std::string> strs;
  std::string text = std::string(buffer);
  boost::split(strs, text , boost::is_any_of("\n"));
  bool pop = false;
  // strs.front().empty() iff \n in front of text
  // strs.back().empty() iff \n in back of text
  for ( auto &s : strs )
  {
    size_t found = s.find_last_of("\r");
    if (found == std::string::npos)
      continue;
    s = s.substr(found+1);
    if(!strs.front().empty())
      pop = true;
  }
  std::string extra = "";
  while (strs.back().empty() && !strs.empty())
  {
    strs.pop_back();
    extra += "\n";
  }
  while (strs.front().empty() && !strs.empty() && strs.front() != strs.back())
  {
    strs.push_back(strs.front());
    strs.pop_front();
  }
  switch (corner)
  {
    case gl::BOTTOM_RIGHT: 
      if (pop) 
        while(!gl::bottomright.empty() && gl::bottomright.front().text.find_first_of("\n") == std::string::npos )
          gl::bottomright.pop_front();
      for ( auto &s : strs )
        gl::bottomright.push_front(gl::TextBlock(extra + s, c));
      break;
    case gl::BOTTOM_LEFT:
      if (pop)
        while(!gl::bottomleft.empty()  && gl::bottomright.front().text.find_first_of("\n") == std::string::npos )
          gl::bottomleft.pop_front();
      for ( auto &s : strs )
        gl::bottomleft.push_front(gl::TextBlock(s + extra, c));
      break;
    case gl::TOP_RIGHT: 
      if (pop)
        while(!gl::topright.empty()    && gl::bottomright.back().text.find_first_of("\n") == std::string::npos )
          gl::topright.pop_back();
      for ( auto &s : strs )
        gl::topright.push_back(gl::TextBlock(extra + s, c));
      break;
    case gl::TOP_LEFT:  // FALL THROUGH
    default: 
      if (pop)
        while(!gl::topleft.empty()     && gl::bottomright.back().text.find_first_of("\n") == std::string::npos )
           gl::topleft.pop_back();
      for ( auto &s : strs )
        gl::topleft.push_back(gl::TextBlock( s + extra, c)); 
      break;
  }
}

// printing routines
#define glprintf_trc(color, fmt, ...) _glprintf(gl::TOP_RIGHT, gl::color, fmt, ##__VA_ARGS__)
#define glprintf_tlc(color, fmt, ...) _glprintf(gl::TOP_LEFT, gl::color, fmt, ##__VA_ARGS__)
#define glprintf_brc(color, fmt, ...) _glprintf(gl::BOTTOM_RIGHT, gl::color, fmt, ##__VA_ARGS__)
#define glprintf_blc(color, fmt, ...) _glprintf(gl::BOTTOM_LEFT, gl::color, fmt, ##__VA_ARGS__)

#define glprintf_tr(fmt, ...) _glprintf(gl::TOP_RIGHT, gl::WHITE, fmt, ##__VA_ARGS__)
#define glprintf_tl(fmt, ...) _glprintf(gl::TOP_LEFT,  gl::WHITE, fmt, ##__VA_ARGS__)
#define glprintf_br(fmt, ...) _glprintf(gl::BOTTOM_RIGHT, gl::WHITE, fmt, ##__VA_ARGS__)
#define glprintf_bl(fmt, ...) _glprintf(gl::BOTTOM_LEFT,  gl::WHITE, fmt, ##__VA_ARGS__)

#endif // GLTEXT_H
