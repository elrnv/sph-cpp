#ifndef SIMWINDOW_H
#define SIMWINDOW_H

#include <QtGui/QOpenGLShaderProgram>
#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLBuffer>
#include "mesh.h"
#include "openglwindow.h"

class SimWindow : public OpenGLWindow
{
public:
    SimWindow();

    void init();
    void render();

private:
    GLuint m_mvpMtxUniform;
    GLuint m_mvInvMtxUniform;
    GLuint m_nmlMtxUniform;
    GLuint m_mvMtxUniform;

    GLuint m_num_vertices;
    GLuint m_num_indices;

		QOpenGLVertexArrayObject *m_vao1;
		QOpenGLBuffer m_posBuffer;
		QOpenGLBuffer m_nmlBuffer;
		QOpenGLBuffer m_idxBuffer;

    QOpenGLShaderProgram *m_program;
    int m_frame;
};

#endif // SIMWINDOW_H
