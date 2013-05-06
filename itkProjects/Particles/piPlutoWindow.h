//
//  piPlutoWindow.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 5/2/13.
//
//

#ifndef __ParticleGuidedRegistration__piPlutoWindow__
#define __ParticleGuidedRegistration__piPlutoWindow__

#include <iostream>

#include <QMainWindow>
#include <QGraphicsScene>
#include <QGraphicsItemGroup>

#include "QGraphicsImageItem.h"
#include "QGraphicsPolygonDrawingInteraction.h"

#include "ui_plutowindow.h"
#include "piPlutoCore.h"

namespace pi {
    typedef QGraphicsImageItem<DataReal> QGraphicsRealImageItem;

    class PlutoWindow: public QMainWindow {
        Q_OBJECT
    public:
        PlutoWindow(QWidget* parent = NULL);
        virtual ~PlutoWindow();
        
        void centerToDesktop();
        void connectSignals();
        
    public slots:
        void on_actionOpen_triggered();
        void on_actionReset_triggered();
        void on_actionStart_triggered();
        void on_actionStep_triggered();
        
        void on_actionLoad_triggered();
        void on_actionSave_triggered();
        
        void flipImages();

    private:
        Ui_PlutoMain _ui;
        QGraphicsScene _scene;
        QGraphicsScene _miniScene;

        RealImage2Vector _images;

        QGraphicsPolygonDrawingInteraction<QGraphicsRealImageItem> _interaction;
        QGraphicsItemGroup* _imageGroup;

        std::vector<QGraphicsRealImageItem*> _imageItems;
        QGraphicsRealImageItem* _patchItem;
        QGraphicsPolygonItem* _test;
        PlutoCore _core;
    };
}

#endif /* defined(__ParticleGuidedRegistration__piPlutoWindow__) */
