TARGET = sim

SOURCES += \
  main.cpp \
  openglwindow.cpp \
  simwindow.cpp \
  scene.cpp \
  mesh.cpp \
  material.cpp \
  primitive.cpp

HEADERS += \
  openglwindow.h \
  simwindow.h \
  scene.h \
  mesh.h \
  material.h \
  primitive.h \
  util.h

QT_CONFIG -= no-pkg-config
CONFIG += link_pkgconfig
PKGCONFIG += eigen3

INCLUDEPATH += /opt/local/include
LIBS += -L/opt/local/lib -lassimp

RESOURCES += resources.qrc

OTHER_FILES += \
  shader.vert \
	shader.frag \
  phong.vert \
  phong.frag \
  normals.geom
