QT += opengl

TARGET = sim

SOURCES += \
  main.cpp \
  gltext.cpp \
  settings.cpp \
  openglwindow.cpp \
  dynamicsmanager.cpp \
  simwindow.cpp \
  mesh.cpp \
  pointcloud.cpp \
  material.cpp \
  uniformbuffer.cpp \
  shadermanager.cpp \
  fluid.cpp \
  boundary.cpp \
  particle.cpp \
  sphgrid.cpp \
  glmesh.cpp \
  glpointcloud.cpp \
  types.cpp

HEADERS += \
  gltext.h \
  settings.h \
  dynparams.h \
  openglwindow.h \
  simwindow.h \
  mesh.h \
  pointcloud.h \
  geometrymanager.h \
  primitive.h \
  material.h \
  materialmanager.h \
  uniformbuffer.h \
  shadermanager.h \
  fluid.h \
  fluiddata.h \
  sphgrid.h \
  boundary.h \
  dynamicsmanager.h \
  glmesh.h \
  glpointcloud.h \
  glprimitive.h \
  particle.h \
  eigen.h \
  extramath.h \
  util.h \
  kernel.h \
  types.h

DESTDIR = ../bin

# directory containing scene configurations
QMAKE_CXXFLAGS += -DCONFIGDIR=/Users/egor/proj/sph/data

# compile with latest c++ specs
CONFIG += c++11

# create a Debug and Release makefiles
CONFIG += debug_and_release

MACPORTS_INCLUDEPATH = /opt/local/include
macx {
  QMAKE_MAC_SDK = macosx10.10
  QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.10
  INCLUDEPATH += $${MACPORTS_INCLUDEPATH}/eigen3
  INCLUDEPATH += $${MACPORTS_INCLUDEPATH}

# Opne Asset Import Library include path
  INCLUDEPATH += $${MACPORTS_INCLUDEPATH}/assimp
# INCLUDEPATH += /Users/egor/proj/assimp/include
# LIBS += -L/Users/egor/proj/assimp/lib

  LIBS += -L/opt/local/lib  # macports libs path
}

LIBS += -lassimp            # 3D asset loading library
LIBS += -ltbb               # concurrency library
LIBS += -lboost_system-mt
LIBS += -lconfig++          # lib to read config files

# LIBS += /opt/intel/lib/libiomp5.a

# MKL setup
#MKLROOT = /opt/intel/mkl
#INCLUDEPATH += $${MKLROOT}/include 
#LIBS += $${MKLROOT}/lib/libmkl_intel_lp64.a $${MKLROOT}/lib/libmkl_core.a $${MKLROOT}/lib/libmkl_intel_thread.a
#LIBS += -L/opt/intel/lib/ -L$${MKLROOT}/lib/ -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread
#LIBS += -ldl -lpthread -lm
#QMAKE_CXXFLAGS += -DMKL_LP64 -m64

# use the following to boost performance
QMAKE_CXXFLAGS += -DBOOST_DISABLE_ASSERTS
#QMAKE_CXXFLAGS += -DREPORT_DENSITY_VARIATION

# show warnings if a function is not defined when inlined
QMAKE_CXXFLAGS += -stdlib=libc++ -Winline -Wno-unused-parameter
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3 -DNDEBUG

RESOURCES += resources.qrc
