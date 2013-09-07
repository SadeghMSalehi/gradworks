//
//  piGroupSimul.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/28/13.
//
//

#ifndef __ParticleGuidedRegistration__piGroupSimul__
#define __ParticleGuidedRegistration__piGroupSimul__

#include <iostream>
#include <QObject>
#include <QVector>
#include "ui_simul2d.h"
#include "piImageDef.h"
#include "QGraphicsImageItem.h"
#include "QGraphicsEventItem.h"
#include "QGraphicsParticleItems.h"
#include <QtConcurrentRun>
#include <QFuture>
#include <QFutureWatcher>
#include <QBasicTimer>

namespace pi {
    class ParticleSystemSolver;
    class ParticleSystem;
    class QGraphicsImageView;
    typedef QGraphicsImageItem<float> QGraphicsRealImageItem;
    typedef QGraphicsImageItem<double> QGraphicsDoubleImageItem;

    class GroupSimul: public QObject, public QGraphicsEllipseEventItem::Listener {
        Q_OBJECT

    public:
        GroupSimul();
        virtual ~GroupSimul();
        
        void setupUi(Ui_Simul2D* ui);
        void setupParticles(ParticleSystemSolver* solver);
        void threadedRun();

    public slots:
        void loadConfig();
        void spreadParticles();
        void computeWarpedImages();
        void computeParticleWarp();
        void computeImageBspline();
        void showWarpedImages();
        void saveTransform();

        void startRun();
        void stopRun();
        
        void showGrid(bool);
        void showParticles(bool);
        void createDistanceMap();
        void showAttributes();
        void showAttribute(int particleId);

        void setResolutionLevel();
        
    protected:
        virtual void timerEvent(QTimerEvent* event);
        virtual void mousePressed(QGraphicsEllipseEventItem* sender, QGraphicsSceneMouseEvent* event);
        virtual void mouseReleased(QGraphicsEllipseEventItem* sender, QGraphicsSceneMouseEvent* event);
        
    private:
        void connectSignals();
        void initImages();
        void initParticles();
        void testTrainer();

    private:
        Ui_Simul2D* _ui;
        ParticleSystemSolver* _solver;

        int _nSubjs;
        int _nParticles;

        QGraphicsScene _scene, _scene2;
        QVector<RealImage::Pointer> _images;
        QVector<RealImage::Pointer> _warpedImages;
        QVector<LabelImage::Pointer> _labels;
        QVector<QGraphicsRealImageItem*> _imageItems;
        QVector<QGraphicsRealImageItem*> _secondImageItems;

        QVector<QGraphicsParticleItems> _particleGroups;
        QVector<QGraphicsParticleItems> _secondParticleGroups;

        QGraphicsRealImageItem* _particleWindow;

        QBasicTimer m_timer;
        QFuture<void> m_future;
        QFutureWatcher<void> m_futureWatcher;
};
}

#endif /* defined(__ParticleGuidedRegistration__piGroupSimul__) */
