TARGET = sim

SOURCES += \
  main.cpp \
  openglwindow.cpp \
  simwindow.cpp \
  scene.cpp \
  mesh.cpp \
  pointcloud.cpp \
  primitive.cpp \
  material.cpp \
  uniformbuffer.cpp \
  shadermanager.cpp \
  dynamics.cpp \
  glmesh.cpp \
  glpointcloud.cpp \

HEADERS += \
  openglwindow.h \
  simwindow.h \
  scene.h \
  mesh.h \
  pointcloud.h \
  primitive.h \
  material.h \
  uniformbuffer.h \
  shadermanager.h \
  dynamics.h \
  glmesh.h \
  glpointcloud.h \
  glprimitive.h \
  eigen.h \
  util.h \
  kernel.h

# compile with latest c++ specs
CONFIG += c++11

QT_CONFIG -= no-pkg-config
CONFIG += link_pkgconfig
PKGCONFIG += eigen3

INCLUDEPATH += /opt/local/include
LIBS += -L/opt/local/lib -lassimp

# use the following to boost performance
# QMAKE_CXXFLAGS += -DBOOST_DISABLE_ASSERTS

# show warnings if a function is not inlined
QMAKE_CXXFLAGS += -Winline

RESOURCES += resources.qrc

OTHER_FILES += \
  phong.vert \
  phong.frag \
  normals.geom \
  normals.vert \
  particle.vert \
  particle.geom \
  particle.frag

