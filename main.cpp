#include <QtGui/QGuiApplication>
#include "simwindow.h"

int main(int argc, char **argv)
{
    QGuiApplication app(argc, argv);
    SimWindow window;
    window.resize(640, 480);
    window.show();

    return app.exec();
}
