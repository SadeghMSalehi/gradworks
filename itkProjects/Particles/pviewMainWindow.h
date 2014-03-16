//
//  pviewMainWindow.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/15/13.
//
//

#ifndef __ParticleGuidedRegistration__pviewMainWindow__
#define __ParticleGuidedRegistration__pviewMainWindow__

#include <iostream>

#include "QMainWindow"
#include "ui_pviewMainWindow.h"

class QContextMenuEvent;
class QFileSystemModel;
class AIRWindow;

class MainWindow: public QMainWindow {
    Q_OBJECT

public:
    MainWindow(QWidget* parent = NULL);
    ~MainWindow();

public slots:
    void on_treeView_clicked(const QModelIndex& index);
    void on_tableView_doubleClicked(const QModelIndex& index);
    void on_tableView_customContextMenuRequested(const QPoint& pos);
    void on_actionLoadMovingImage_triggered();

private:
    Ui::MainWindow ui;
    QFileSystemModel* m_dirModel;
    QFileSystemModel* m_fileModel;
    AIRWindow* m_AIRWindow;
};

#endif /* defined(__ParticleGuidedRegistration__pviewMainWindow__) */
