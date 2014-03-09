//
//  airSuperZoom.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/30/13.
//
//

#include "airSuperZoom.h"
#include <QGLWidget>
#include <QGraphicsRectItem>
#include <QList>

using namespace std;


namespace air {
    SuperZoom::SuperZoom(QWidget* p): QDialog(p) {
        ui.setupUi(this);

        ui.graphicsView->setScene(&_scene);
        ui.graphicsView->setViewport(new QGLWidget(this));
        ui.graphicsView->setDragMode(QGraphicsView::ScrollHandDrag);

        QList<QAction*> menuActions;
        menuActions.push_back(ui.actionShowSomething);
        ui.toolButton->addActions(menuActions);
        ui.toolButton->setDefaultAction(ui.actionShowSomething);

        connect(ui.actionShowSomething, SIGNAL(triggered()), this, SLOT(actionTriggered()));
    }

    SuperZoom::~SuperZoom() {
        
    }

    void SuperZoom::actionTriggered() {
        QAction* action = dynamic_cast<QAction*>(sender());
        if (action == ui.actionShowSomething) {
            QGraphicsRectItem* item = _scene.addRect(0,0,100,100,QPen(Qt::yellow,3));
            item->setFlag(QGraphicsItem::ItemIsMovable);
            ui.graphicsView->centerOn(item);
            //ui.graphicsView->fitInView(item);
        }
    }

    void SuperZoom::showEvent(QShowEvent *event) {
        
    }

    void SuperZoom::closeEvent(QCloseEvent* event) {
    }
}