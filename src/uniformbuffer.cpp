#include "uniformbuffer.h"

#include <QOpenGLContext>
#include <QtGui/QOpenGLFunctions>

#include <QDebug>

QVector<GLuint> UniformBuffer::m_bindingIndices = QVector<GLuint>();

GLuint UniformBuffer::bindingIndex() 
{
  for (GLuint i = 0, size = m_bindingIndices.size(); i < size; i++) 
  {
    if (m_bindingIndices.at(i) != i) 
    {
      m_bindingIndices.insert(i, i);
      return i;
    }
  }
  m_bindingIndices.append( m_bindingIndices.size() );
  return m_bindingIndices.size();
}

UniformBuffer::UniformBuffer()
 : m_bufferId(0)
 , m_bindingIdx(-1)
 , m_glContext(NULL)
 , m_inited(false)
 //, glBindBuffer(NULL)
 , glBindBufferBase(NULL)
 , glBindBufferRange(NULL)
 //, glBufferData(NULL)
 //, glDeleteBuffers(NULL)
 //, glGenBuffers(NULL)
 , glGetActiveUniformBlockiv(NULL), glGetActiveUniformsiv(NULL)
 , glGetUniformBlockIndex(NULL), glGetUniformIndices(NULL)
 , glUniformBlockBinding(NULL)
{
}


UniformBuffer::UsagePattern UniformBuffer::usagePattern() const
{
  return m_usage_pattern;
}

void UniformBuffer::setUsagePattern(UniformBuffer::UsagePattern value)
{
  m_usage_pattern = value;
}

bool UniformBuffer::create() 
{
  const QOpenGLContext* m_glContext = QOpenGLContext::currentContext();

  if (m_glContext) 
  {
    if (!m_inited && !initFunctionPointers(m_glContext))
    {
      qDebug("Cannot find Uniform Buffer Objects related functions");
      return false;
    }

    GLuint tmpBufferId = 0;
    m_glContext->functions()->glGenBuffers(1, &tmpBufferId);

    if (tmpBufferId) 
    {
      m_bufferId = tmpBufferId;
      this->m_glContext = m_glContext;

      return true;
    }
    else
    {
      qDebug("Invalid buffer Id");
    }
  }

  qDebug("Could not retrieve buffer");

  return false;
}

bool UniformBuffer::isCreated() const 
{
  return (bool)(m_bufferId != 0);
}

void UniformBuffer::destroy() 
{
  if ( m_bufferId || m_glContext ) 
  {
    m_glContext->functions()->glDeleteBuffers(1, &m_bufferId);
  }
#ifndef QT_NO_DEBUG
  else { qWarning("UniformBuffer::destroy(): buffer already destroyed."); }
#endif
  m_bufferId = 0;
  m_glContext = NULL;
}

bool UniformBuffer::bind() 
{
  if (!isCreated()) 
  {
    qDebug("Buffer not created");
    return false;
  }

  m_glContext->functions()->glBindBuffer(GL_UNIFORM_BUFFER, m_bufferId);

  return true;
}

int UniformBuffer::bindToIndex() 
{
  m_bindingIdx = bindingIndex();
  glBindBufferBase(GL_UNIFORM_BUFFER, m_bindingIdx, m_bufferId);
  if (glGetError() == GL_INVALID_VALUE)
  {
    qDebug() << "OpenGL ERROR: GL_INVALID_VALUE, when binding uniform buffer to index.";
    return -1;
  }

  return m_bindingIdx;
}

bool UniformBuffer::releaseFromIndex() 
{
  glBindBufferBase(GL_UNIFORM_BUFFER, m_bindingIdx, 0);
  if (glGetError() == GL_INVALID_VALUE)
  {
    qDebug() << "OpenGL ERROR: GL_INVALID_VALUE, when unbinding uniform buffer from index.";
    return false;
  }
  m_bindingIdx = -1;

  return true;
}

bool UniformBuffer::bindToProg(GLuint progId, QString uniformName) 
{
  GLuint tmpBlockIdx = glGetUniformBlockIndex(progId, uniformName.toUtf8());

  if (tmpBlockIdx == GL_INVALID_INDEX) 
  {
    qDebug() << QString("Could not find block index of block named: %1").arg(uniformName);

    return false;
  }

  //GLint tmpBlockSize;
  //glGetActiveUniformBlockiv(progId, tmpBlockIdx, GL_UNIFORM_BLOCK_DATA_SIZE, &tmpBlockSize);

  glUniformBlockBinding(progId, tmpBlockIdx, m_bindingIdx);

  if (glGetError() == GL_INVALID_VALUE)
  {
    qDebug() << "OpenGL ERROR: GL_INVALID_VALUE, when binding uniform buffer to shader program.";
    return false;
  }

  //UBOInfo info;
  //info.progId = progId;
  //info.uniformName = uniformName;
  //info.blockIdx = tmpBlockIdx;
  //info.blockSize = tmpBlockSize;
  //info.bindingIdx = tmpBindingIdx;

  //m_UBOInfos.append(info);

  return true;
}

void UniformBuffer::release() 
{
  m_glContext->functions()->glBindBuffer(GL_UNIFORM_BUFFER, 0);
}

GLuint UniformBuffer::bufferId() const 
{
  return m_bufferId;
}


// returns offset + count if data is written
int UniformBuffer::write(int offset, const void *data, int count) 
{
  if (!isCreated())
  {
#ifndef QT_NO_DEBUG
    qWarning("UnformBuffer::write(): buffer not created");
#endif
    return offset;
  }

  m_glContext->functions()->glBufferSubData(GL_UNIFORM_BUFFER, offset, count, data);
  return offset + count;
}

void UniformBuffer::allocate(const void *data, int count) 
{
  if (!isCreated())
  {
#ifndef QT_NO_DEBUG
    qWarning("UnformBuffer::allocate(): buffer not created");
#endif
    return;
  }

  m_glContext->functions()->glBufferData(GL_UNIFORM_BUFFER, count, data, m_usage_pattern);
}

bool UniformBuffer::initFunctionPointers(const QOpenGLContext* m_glContext)
{
  //glBindBuffer = (PFNGLBINDBUFFERPROC)m_glContext->getProcAddress("glBindBuffer");
  glBindBufferBase = (PFNGLBINDBUFFERBASEPROC)m_glContext->getProcAddress("glBindBufferBase");
  glBindBufferRange = (PFNGLBINDBUFFERRANGEPROC)m_glContext->getProcAddress("glBindBufferRange");
  //glBufferData = (PFNGLBUFFERDATAPROC)m_glContext->getProcAddress("glBufferData");
  //glBufferSubData = (PFNGLBUFFERSUBDATAPROC)m_glContext->getProcAddress("glBufferSubData");
  //glDeleteBuffers = (PFNGLDELETEBUFFERSPROC)m_glContext->getProcAddress("glDeleteBuffers");
  //glGenBuffers = (PFNGLGENBUFFERSPROC)m_glContext->getProcAddress("glGenBuffers");
  glGetActiveUniformBlockiv = (PFNGLGETACTIVEUNIFORMBLOCKIVPROC)m_glContext->getProcAddress("glGetActiveUniformBlockiv");
  glGetActiveUniformsiv = (PFNGLGETACTIVEUNIFORMSIVPROC)m_glContext->getProcAddress("glGetActiveUniformsiv");
  glGetUniformBlockIndex = (PFNGLGETUNIFORMBLOCKINDEXPROC)m_glContext->getProcAddress("glGetUniformBlockIndex");
  glGetUniformIndices = (PFNGLGETUNIFORMINDICESPROC)m_glContext->getProcAddress("glGetUniformIndices");
  glUniformBlockBinding = (PFNGLUNIFORMBLOCKBINDINGPROC)m_glContext->getProcAddress("glUniformBlockBinding");

  if(//!glBindBuffer ||
      !glBindBufferBase ||
      !glBindBufferRange ||
      //!glBufferData ||
      //!glBufferSubData ||
      //!glDeleteBuffers ||
      //!glGenBuffers ||
      !glGetActiveUniformBlockiv ||
      !glGetActiveUniformsiv ||
      !glGetUniformBlockIndex ||
      !glGetUniformIndices ||
      !glUniformBlockBinding)
  {
    qDebug("Could not init function pointers");
    return false;
  }

  return true;
}

