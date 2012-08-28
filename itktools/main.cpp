#include "itkImageCommon.h"
#include "mainwindow.h"
#include <QtGui>
#include <QApplication>

using namespace std;
using namespace itk;

typedef Image<short, 3> ImageType;

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
}
