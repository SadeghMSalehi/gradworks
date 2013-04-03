//
//  dualImageViewer.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/3/13.
//
//

#ifndef __ParticleGuidedRegistration__dualImageViewer__
#define __ParticleGuidedRegistration__dualImageViewer__

#include <iostream>

#include <QMainWindow>
#include <QGraphicsScene>

#include "ui_dualImageViewer.h"

class DualImageViewer: public QMainWindow {
    Q_OBJECT
public:
    DualImageViewer(QWidget* parent = NULL, Qt::WindowFlags = 0);
    virtual ~DualImageViewer();

public slots:
    void viewMouseEvent(QMouseEvent* mouse);
    void clearScene();
    
private:
    Ui::DualViewWindow ui;

    void setupUi();
    void connectSignals();
};

#endif /* defined(__ParticleGuidedRegistration__dualImageViewer__) */
