//
//  pviewMainWindow.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/15/13.
//
//

#include "pviewMainWindow.h"
#include "QFileSystemModel"
#include "pviewImageViewer.h"
#include "QContextMenuEvent"

MainWindow::MainWindow(QWidget* parent) {
    m_imageViewer = NULL;

    ui.setupUi(this);
    m_dirModel = new QFileSystemModel();
    m_fileModel = new QFileSystemModel();
    m_dirModel->setFilter(QDir::AllDirs | QDir::NoDotAndDotDot);
    m_dirModel->sort(0);
    ui.treeView->setModel(m_dirModel);
    m_dirModel->setRootPath("/NIRAL/work/joohwi");
    ui.treeView->setRootIndex(m_dirModel->index("/NIRAL/work/"));

    ui.treeView->hideColumn(1);
    ui.treeView->hideColumn(2);
    ui.treeView->hideColumn(3);
    ui.treeView->setIndentation(20);
    ui.treeView->setSortingEnabled(true);
    ui.treeView->setWindowTitle("Directories");

    m_fileModel->setRootPath("/NIRAL/work/joohwi");
    ui.tableView->setModel(m_fileModel);
    ui.tableView->setRootIndex(m_fileModel->index("/NIRAL/work/joohwi"));
    ui.tableView->verticalHeader()->hide();
    ui.tableView->setSortingEnabled(true);
    ui.tableView->setSelectionBehavior(QAbstractItemView::SelectRows);
    ui.tableView->resizeColumnToContents(0);
    ui.tableView->setContextMenuPolicy(Qt::CustomContextMenu);
}

MainWindow::~MainWindow() {

}

void MainWindow::on_treeView_clicked(const QModelIndex& index) {
    QString sPath = m_dirModel->fileInfo(index).absoluteFilePath();
    ui.tableView->setRootIndex(m_fileModel->setRootPath(sPath));
}

void MainWindow::on_tableView_doubleClicked(const QModelIndex &index) {
    QFileInfo fileInfo = m_fileModel->fileInfo(index);
    QString sPath = fileInfo.absoluteFilePath();
    if (fileInfo.isDir()) {
        ui.tableView->setRootIndex(m_fileModel->index(sPath));
        ui.treeView->setCurrentIndex(m_dirModel->index(sPath));
    } else {
        if (m_imageViewer == NULL) {
            m_imageViewer = new ImageViewer();
        }
        if (m_imageViewer->isHidden()) {
            m_imageViewer->LoadImage(sPath);
            m_imageViewer->show();
        } else {
            m_imageViewer->LoadImage(sPath);
        }
    }
}

void MainWindow::on_tableView_customContextMenuRequested(const QPoint& pos) {
    if (m_imageViewer != NULL && m_imageViewer->isVisible()) {
        QMenu menu(this);
        menu.addAction(ui.actionLoadMovingImage);
        menu.exec(ui.tableView->mapToGlobal(pos));

        QModelIndex cellIdx = ui.tableView->indexAt(pos);
        ui.tableView->setCurrentIndex(cellIdx);
    }
}

void MainWindow::on_actionLoadMovingImage_triggered() {
    QModelIndex idx = ui.tableView->currentIndex();
    QFileInfo fileInfo = m_fileModel->fileInfo(idx);
    m_imageViewer->LoadImage(fileInfo.absoluteFilePath());
}