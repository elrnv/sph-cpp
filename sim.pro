QT += opengl

TARGET = sim

SOURCES += \
  main.cpp \
  gltext.cpp \
  settings.cpp \
  openglwindow.cpp \
  simwindow.cpp \
  scene.cpp \
  mesh.cpp \
  pointcloud.cpp \
  primitive.cpp \
  material.cpp \
  uniformbuffer.cpp \
  shadermanager.cpp \
#  quantityprocessor.cpp \
  fluid.cpp \
  dynamics.cpp \
  glmesh.cpp \
  glpointcloud.cpp

HEADERS += \
  gltext.h \
  settings.h \
  dynparams.h \
  openglwindow.h \
  simwindow.h \
  scene.h \
  mesh.h \
  pointcloud.h \
  primitive.h \
  material.h \
  uniformbuffer.h \
  shadermanager.h \
  fluidmanager.h \
#  quantityprocessor.h \
  fluid.h \
  dynamics.h \
  glmesh.h \
  glpointcloud.h \
  glprimitive.h \
  particle.h \
  eigen.h \
  util.h \
  kernel.h

DESTDIR = ../bin

# compile with latest c++ specs
CONFIG += c++11

INCLUDEPATH += /opt/local/include/eigen3
INCLUDEPATH += /opt/local/include
LIBS += -L/opt/local/lib 
LIBS += -lassimp
LIBS += -lboost_system-mt

# LIBS += /opt/intel/lib/libiomp5.a

# link to libconfig library to read config files
LIBS += -lconfig++

# MKL setup
#MKLROOT = /opt/intel/mkl
#INCLUDEPATH += $${MKLROOT}/include 
#LIBS += $${MKLROOT}/lib/libmkl_intel_lp64.a $${MKLROOT}/lib/libmkl_core.a $${MKLROOT}/lib/libmkl_intel_thread.a
#LIBS += -L/opt/intel/lib/ -L$${MKLROOT}/lib/ -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread
#LIBS += -ldl -lpthread -lm
#QMAKE_CXXFLAGS += -DMKL_LP64 -m64

# use the following to boost performance
QMAKE_CXXFLAGS += -DBOOST_DISABLE_ASSERTS

# show warnings if a function is not inlined
QMAKE_CXXFLAGS += -stdlib=libc++ -Winline -Wno-unused-parameter

# directory containing scene configurations
QMAKE_CXXFLAGS += -DCONFIGDIR=/Users/egor/proj/sph/data

RESOURCES += resources.qrc
