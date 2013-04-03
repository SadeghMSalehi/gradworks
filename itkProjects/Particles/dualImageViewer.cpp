//
//  dualImageViewer.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/3/13.
//
//

#include "dualImageViewer.h"
#include <QMouseEvent>

using namespace std;

DualImageViewer::DualImageViewer(QWidget* parent, Qt::WindowFlags flags): QMainWindow(parent, flags) {
    ui.setupUi(this);
    setupUi();
    connectSignals();
}

DualImageViewer::~DualImageViewer() {
    
}

void DualImageViewer::setupUi() {
    QGraphicsScene* scene = new QGraphicsScene(this);
    ui.view1->setScene(scene);
    ui.view2->setScene(scene);

    scene->setSceneRect(ui.view1->viewport()->contentsRect());
}

void DualImageViewer::connectSignals() {
    connect(ui.actionClear, SIGNAL(triggered()), SLOT(clearScene()));
    connect(ui.view1, SIGNAL(mousePressed(QMouseEvent*)), SLOT(viewMouseEvent(QMouseEvent*)));
    connect(ui.view1, SIGNAL(mouseMoved(QMouseEvent*)), SLOT(viewMouseEvent(QMouseEvent*)));
    connect(ui.view1, SIGNAL(mouseReleased(QMouseEvent*)), SLOT(viewMouseEvent(QMouseEvent*)));
}

void DualImageViewer::clearScene() {
    ui.view1->scene()->clear();
}

void DualImageViewer::viewMouseEvent(QMouseEvent *e) {
    QGraphicsScene* scene = ui.view1->scene();
    QPointF pos = ui.view1->mapToScene(e->pos());
    scene->addEllipse(pos.x()-2.5, pos.y()-2.5, 5, 5, QPen(Qt::yellow), QBrush(Qt::yellow));
}