#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QMessageBox>
#include <QProcess>
#include <QFileDialog>

#include <iostream>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    connect(&_script, SIGNAL(readyReadStandardOutput()), this, SLOT(readScriptOutput()));
    connect(&_script, SIGNAL(stateChanged(QProcess::ProcessState)), this, SLOT(scriptStateChanged(QProcess::ProcessState)));
}

MainWindow::~MainWindow()
{
    delete ui;
}


std::string MainWindow::findProgram(std::string name) {
    using namespace std;

    string command = "which " + name;

    QProcess finder(this);
    finder.setProcessChannelMode(QProcess::MergedChannels);
    finder.start("/bin/sh", QStringList() << "-c" << QString(command.c_str()));

    if (!finder.waitForStarted()) {
        return "";
    }

    QString pathToName;
    if (finder.waitForReadyRead()) {
        pathToName = finder.readAll().constData();
        pathToName = pathToName.simplified();
    }

    finder.waitForFinished();
    return pathToName.toStdString();
}


void MainWindow::runScript() {
    using namespace std;

    if (_script.state() != QProcess::NotRunning) {
        return;
    }

    cout << "'" << findProgram("python") << "'" << endl;
    cout << "'" << findProgram("ImageMath") << "'" << endl;


//    _script.start("/usr/bin/python", QStringList() << "/tools/gradworks/test/qt/ProcessRunner/test.py");
}

void MainWindow::stopScript() {
    if (_script.state() != QProcess::Running) {
        return;
    }
    _script.kill();
}

void MainWindow::readScriptOutput() {
    using namespace std;

    cout << _script.readAll().constData() << flush;
}

void MainWindow::scriptStateChanged(QProcess::ProcessState state) {
    using namespace std;

    if (state == QProcess::NotRunning) {
        cout << "script has finished ..." << endl;
    } else if (state == QProcess::Running) {
        cout << "script has begun ..." << endl;
    }
}


void MainWindow::addRow() {
    int itemRow;
    itemRow = ui->tableWidget->row(ui->tableWidget->selectedItems().at(0));
    if (itemRow < 0) {
        itemRow = ui->tableWidget->rowCount();
        ui->tableWidget->QTableWidget::insertRow(itemRow);
    } else {
        ui->tableWidget->QTableWidget::insertRow(itemRow+1);
    }
}

void MainWindow::deleteCurrentRow() {
    int itemRow = ui->tableWidget->row(ui->tableWidget->selectedItems().at(0));
    if (itemRow < 0) {
        itemRow = ui->tableWidget->rowCount();
        ui->tableWidget->QTableWidget::removeRow(itemRow-1);
    } else {
        ui->tableWidget->QTableWidget::removeRow(itemRow);
    }
}

void MainWindow::editCell(int row, int column) {
    if (column == 1) {
        QString item = QFileDialog::getOpenFileName(this, "Open image", QString(), "NRRD Images (*.nrrd *.nhdr)");
        QTableWidgetItem *itemCell = new QTableWidgetItem(tr("%1").arg(item));
        ui->tableWidget->setItem(row, column, itemCell);
    }
}
