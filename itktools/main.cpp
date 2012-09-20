#include "itkMyCore.h"
#include "mainwindow.h"
#include <QtGui>
#include <QApplication>

using namespace std;
using namespace itk;

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
}
