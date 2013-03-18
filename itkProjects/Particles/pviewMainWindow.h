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

class QFileSystemModel;

class MainWindow: public QMainWindow {
    Q_OBJECT

public:
    MainWindow(QWidget* parent = NULL);
    ~MainWindow();

public slots:
    void on_treeView_clicked(const QModelIndex& index);
    void on_tableView_doubleClicked(const QModelIndex& index);

private:
    Ui::MainWindow ui;
    QFileSystemModel* m_dirModel;
    QFileSystemModel* m_fileModel;
};

#endif /* defined(__ParticleGuidedRegistration__pviewMainWindow__) */
