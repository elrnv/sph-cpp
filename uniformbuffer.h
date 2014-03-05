#ifndef UNIFORMBUFFER_H
#define UNIFORMBUFFER_H

#include <QtOpenGL/qgl.h>

struct UBOInfo
{
  QString uniformName;
  GLuint progId;
  GLuint blockIdx;
  GLint blockSize;
  GLuint bindingIdx;
};


class UniformBuffer
{
public:
  UniformBuffer();

  enum UsagePattern
  {
    StreamDraw   = 0x88E0, // GL_STREAM_DRAW
    StreamRead   = 0x88E1, // GL_STREAM_READ
    StreamCopy   = 0x88E2, // GL_STREAM_COPY
    StaticDraw   = 0x88E4, // GL_STATIC_DRAW
    StaticRead   = 0x88E5, // GL_STATIC_READ
    StaticCopy   = 0x88E6, // GL_STATIC_COPY
    DynamicDraw  = 0x88E8, // GL_DYNAMIC_DRAW
    DynamicRead  = 0x88E9, // GL_DYNAMIC_READ
    DynamicCopy  = 0x88EA  // GL_DYNAMIC_COPY
  };

  UniformBuffer::UsagePattern usagePattern() const;
  void setUsagePattern(UniformBuffer::UsagePattern value);

  bool create();
  bool isCreated() const;

  void destroy();

  bool bind();
  int bindToIndex(); // return found m_bindingIdx
  bool releaseFromIndex();
  bool bindToProg(GLuint progId, QString uniformName);
  void release();

  GLuint bufferId() const;

  //    bool read(int offset, void *data, int count);
  int write(int offset, const void *data, int count);

  void allocate(const void *data, int count);
  inline void allocate(int count) { allocate(0, count); }


protected:
  static GLuint bindingIndex();
  static QVector<GLuint> m_bindingIndices;

  bool initFunctionPointers(const QOpenGLContext* m_glContext);

  GLuint m_bufferId;
  GLuint m_bindingIdx;

  const QOpenGLContext* m_glContext;

  bool m_inited;

  UniformBuffer::UsagePattern m_usage_pattern;

  QVector<UBOInfo> m_UBOInfos;

  //PFNGLBINDBUFFERPROC glBindBuffer;
  PFNGLBINDBUFFERBASEPROC glBindBufferBase;
  PFNGLBINDBUFFERRANGEPROC glBindBufferRange;
  //PFNGLBUFFERDATAPROC glBufferData;
  //PFNGLBUFFERSUBDATAPROC glBufferSubData;
  //PFNGLDELETEBUFFERSPROC glDeleteBuffers;
  //PFNGLGENBUFFERSPROC glGenBuffers;
  PFNGLGETACTIVEUNIFORMBLOCKIVPROC glGetActiveUniformBlockiv;
  PFNGLGETACTIVEUNIFORMSIVPROC glGetActiveUniformsiv;
  PFNGLGETUNIFORMBLOCKINDEXPROC glGetUniformBlockIndex;
  PFNGLGETUNIFORMINDICESPROC glGetUniformIndices;
  PFNGLUNIFORMBLOCKBINDINGPROC glUniformBlockBinding;
};

#endif // UNIFORMBUFFER_H
