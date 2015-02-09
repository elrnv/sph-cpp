QT += opengl

TARGET = sph

SRCDIR = src
INCLUDEDIR = src

SOURCES += \
  $${SRCDIR}/main.cpp \
  $${SRCDIR}/gltext.cpp \
  $${SRCDIR}/settings.cpp \
  $${SRCDIR}/openglwindow.cpp \
  $${SRCDIR}/sphgrid.cpp \
  $${SRCDIR}/dynamicsmanager.cpp \
  $${SRCDIR}/simwindow.cpp \
  $${SRCDIR}/mesh.cpp \
  $${SRCDIR}/pointcloud.cpp \
  $${SRCDIR}/material.cpp \
  $${SRCDIR}/uniformbuffer.cpp \
  $${SRCDIR}/shadermanager.cpp \
  $${SRCDIR}/fluid.cpp \
  $${SRCDIR}/boundary.cpp \
  $${SRCDIR}/particle.cpp \
  $${SRCDIR}/glmesh.cpp \
  $${SRCDIR}/glpointcloud.cpp \
  $${SRCDIR}/types.cpp

HEADERS += \
  $${INCLUDEDIR}/gltext.h \
  $${INCLUDEDIR}/settings.h \
  $${INCLUDEDIR}/dynparams.h \
  $${INCLUDEDIR}/openglwindow.h \
  $${INCLUDEDIR}/simwindow.h \
  $${INCLUDEDIR}/mesh.h \
  $${INCLUDEDIR}/pointcloud.h \
  $${INCLUDEDIR}/geometrymanager.h \
  $${INCLUDEDIR}/primitive.h \
  $${INCLUDEDIR}/material.h \
  $${INCLUDEDIR}/materialmanager.h \
  $${INCLUDEDIR}/uniformbuffer.h \
  $${INCLUDEDIR}/shadermanager.h \
  $${INCLUDEDIR}/fluid.h \
  $${INCLUDEDIR}/fluiddata.h \
  $${INCLUDEDIR}/sphgrid.h \
  $${INCLUDEDIR}/boundary.h \
  $${INCLUDEDIR}/dynamicsmanager.h \
  $${INCLUDEDIR}/glmesh.h \
  $${INCLUDEDIR}/glpointcloud.h \
  $${INCLUDEDIR}/glprimitive.h \
  $${INCLUDEDIR}/particle.h \
  $${INCLUDEDIR}/eigen.h \
  $${INCLUDEDIR}/extramath.h \
  $${INCLUDEDIR}/util.h \
  $${INCLUDEDIR}/kernel.h \
  $${INCLUDEDIR}/types.h

RESOURCES += resources.qrc

# directory containing scene configurations
QMAKE_CXXFLAGS += -DCONFIGDIR=/Users/egor/proj/sph/data

# compile with latest c++ specs
CONFIG += c++11

# create a Debug and Release makefiles
CONFIG += debug_and_release

CONFIG(release, debug|release) {
  DESTDIR = build/release
} else {
  DESTDIR = build/debug
}

OBJECTS_DIR = $${DESTDIR}/.obj
MOC_DIR =     $${DESTDIR}/.moc
RCC_DIR =     $${DESTDIR}/.qrc
UI_DIR =      $${DESTDIR}/.ui

MACPORTS_INCLUDEPATH = /opt/local/include
HOMEBREW_INCLUDEPATH = /usr/local/include

macx {
  QMAKE_MAC_SDK = macosx10.10
  QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.10
  INCLUDEPATH += $${HOMEBREW_INCLUDEPATH}/eigen3
  INCLUDEPATH += $${HOMEBREW_INCLUDEPATH}

# Opne Asset Import Library include path
#  INCLUDEPATH += $${HOMEBREW_INCLUDEPATH}/assimp

# INCLUDEPATH += /Users/egor/proj/assimp/include
# LIBS += -L/Users/egor/proj/assimp/lib

  LIBS += -L/opt/local/lib  # macports libs path
  LIBS += -L/usr/local/lib  # homebrew libs path
}

# this is needed so eigen can find our eigen plugins
INCLUDEPATH += $${SRCDIR}

LIBS += -lassimp            # 3D asset loading library
#LIBS += -ltbb               # concurrency library
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
#QMAKE_CXXFLAGS += -DREPORT_DENSITY_VARIATION

QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3 -DNDEBUG -DBOOST_DISABLE_ASSERTS
