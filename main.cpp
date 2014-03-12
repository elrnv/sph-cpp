#include <QtGui/QGuiApplication>
#include "simwindow.h"

int main(int argc, char **argv)
{
  QGuiApplication app(argc, argv);
  SimWindow window;
  window.resize(800, 600);
  window.show();

  return app.exec();
}
