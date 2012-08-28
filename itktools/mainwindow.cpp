#include "mainwindow.h"
#include "QMessageBox"

MainWindow::MainWindow(QWidget *parent) :
		QMainWindow(parent) {
	ui.setupUi(this);
//	QObject::connect(ui.action_New, SIGNAL(triggered(bool)), this,
//			SLOT(newMenu(void)));
//	QObject::connect(ui.action_Close, SIGNAL(triggered(bool)), this,
//			SLOT(quitMenu(void)));
	if (!connect(this->ui.action_New, SIGNAL(triggered(bool)), this,
			SLOT(newMenu(void)))) {
		newMenu();
	}
}

MainWindow::~MainWindow() {

}

void MainWindow::newMenu() {
	QMessageBox::information(this, tr("New?"), tr("Do you want to new?"),
			QMessageBox::Ok, QMessageBox::Ok);
}

void MainWindow::quitMenu() {
	QMessageBox::information(this, tr("Quit?"), tr("Do you want to quit?"),
			QMessageBox::Ok, QMessageBox::Ok);
}
