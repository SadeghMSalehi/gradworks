#include "mainwindow.h"
#include "QMessageBox"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
	ui.setupUi(this);
	connect(ui.action_New, SIGNAL(triggered(bool)), this, SLOT(sayHello(void)));
}

MainWindow::~MainWindow()
{

}

void MainWindow::sayHello(void) {
	QMessageBox::information(this, tr("Apps"), tr("Hello~"), QMessageBox::Ok, QMessageBox::Ok);
}
