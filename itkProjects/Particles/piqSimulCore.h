//
//  piqSimulCore.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/7/13.
//
//

#ifndef __ParticleGuidedRegistration__piqSimulCore__
#define __ParticleGuidedRegistration__piqSimulCore__

#include <iostream>

#include <QObject>
#include <QVector>

#include "ui_simul2d.h"
#include "piImageDef.h"
#include "piParticle.h"

template <class T> class QGraphicsImageItem;

class QGraphicsScene;
class QGraphicsPixmapItem;
class QGraphicsEllipseItem;
class QGraphicsRectWidget;

namespace pi {
    class ParticleSystem;
    class ParticleSystemSolver;
}

/**
 * This class is responsible for running and executing registration algorithm.
 *
 */
namespace piq {
    class SimulCore: public QObject {
        Q_OBJECT

    public:
        typedef QGraphicsImageItem<pi::RealImage> QRealImageItem;

        SimulCore(QWidget* parent = NULL);
        virtual ~SimulCore();
        void setParticleSolver(pi::ParticleSystemSolver* solver);
        void showAuxImage(int, pi::RealImage::Pointer image);

    public slots:
        void setUi(Ui_Simul2D* ui);
        void openImage(int, QString);
        void openLabel(int, QString);
        void hideParticles(bool show);
        
        void startTrackingMode();
        void alignTrackingTarget();
        void trackingWidgetMoved(QPointF);
        
        void labelOpacityChanged(int value);
        void updateParticles();
        void updateParticles(int n, pi::ParticleVector& particles);

        void run();


        QGraphicsItem* getImageItem(int n);
        QGraphicsScene* scene(int n);

    private:
        void connectSignals();
        QRealImageItem* showImage(int n, pi::RealImage::Pointer image);
        QGraphicsPixmapItem* showLabel(int n, pi::LabelImage::Pointer image);
        void createParticleItems(int id, int n);

    private:
        Ui_Simul2D* _ui;
        QWidget* _parent;

        pi::RealImage::Pointer _image[2];
        pi::LabelImage::Pointer _label[2];

        QGraphicsScene* _scene[2];
        QGraphicsScene* _miniScene;
        
        QGraphicsImageItem<pi::RealImage>* _imageItem[2];
        QGraphicsImageItem<pi::RealImage>* _auxImageItem[2];
        
        QGraphicsPixmapItem* _labelItem[2];
        QVector<QGraphicsEllipseItem*> _particleItem[2];
        QGraphicsRectWidget* _rectItem[2];

        pi::ParticleSystemSolver* _solver;
    };
}
#endif /* defined(__ParticleGuidedRegistration__piqSimulCore__) */
