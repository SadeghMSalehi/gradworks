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

#include "ui_simul2d.h"
#include "piImageDef.h"
#include "piParticle.h"

template <class T> class QGraphicsImageItem;
class QGraphicsScene;
class QGraphicsPixmapItem;
class QParticlesGraphicsItem;

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

    public slots:
        void setUi(Ui_Simul2D* ui);
        void openImage1(QString);
        void openImage2(QString);
        void openLabel1(QString);
        void openLabel2(QString);

        void labelOpacityChanged(int value);
        void updateParticles();
        void updateParticles(int n, pi::ParticleVector& particles);

        void run();


        QGraphicsItem* getImageItem(int n);
        QGraphicsScene* scene(int n);

    private:
        void connectSignals();
        QRealImageItem* showImage(QGraphicsScene* scene, QRealImageItem* item, pi::RealImage::Pointer image);
        QGraphicsPixmapItem* showLabel(QGraphicsScene* scene, QGraphicsPixmapItem* item, pi::LabelImage::Pointer image, QGraphicsItem* parent);

    private:
        Ui_Simul2D* _ui;
        QWidget* _parent;

        pi::RealImage::Pointer _image[2];
        pi::LabelImage::Pointer _label[2];

        QGraphicsScene* _scene[2];
        
        QGraphicsImageItem<pi::RealImage>* _imageItem[2];
        QGraphicsPixmapItem* _labelItem[2];
        QParticlesGraphicsItem* _particles[2];

        pi::ParticleSystemSolver* _solver;
    };
}
#endif /* defined(__ParticleGuidedRegistration__piqSimulCore__) */
