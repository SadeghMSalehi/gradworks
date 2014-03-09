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
#include "piPatchTracking.h"
#include "QGraphicsEventItem.h"
#include "QGraphicsImageItem.h"

template <class T> class QGraphicsImageItem;
typedef QGraphicsImageItem<pi::RealImage::PixelType> QGraphicsRealImageItem;

class QGraphicsScene;
class QGraphicsPixmapItem;
class QGraphicsEllipseItem;
class QGraphicsRectWidget;
class QGraphicsGridItem;
class QAbstractGraphicsShapeItem;
class SimulCore;

namespace pi {
    class ParticleSystem;
    class ParticleSystemSolver;
    class SimulCore;
}

typedef QGraphicsEventItem<QGraphicsEllipseItem> EllipseItem;

/**
 * This class is responsible for running and executing registration algorithm.
 *
 */
namespace piq {
    class ImagePatchInteraction;
    
    class SimulCore: public QObject, public EllipseItem::Listener, public QGraphicsRealImageItem::InteractionType {
        Q_OBJECT

    public:
        enum { ParticleId, ImageId };
        enum SelectionMode { None, Selected };
        enum ImageInteractionMode { NoInteraction, InspectionMode, TrackingMode, SIFTMode };

        SimulCore(QWidget* parent = NULL);
        virtual ~SimulCore();
        void setParticleSolver(pi::ParticleSystemSolver* solver);
        void showAuxImage(int, pi::RealImage::Pointer image, bool autoRange = false);
        void showCompositeImage(pi::RealImage::Pointer, pi::RealImage::Pointer, pi::RealImage::Pointer);
        void showEntropyMeasurementImage();

        pi::RealImage::Pointer getImage(int i);
        QGraphicsRealImageItem* getImageItem(int n);

    public slots:
        void setUi(Ui_Simul2D* ui);
        void openImage(int, QString);
        void openLabel(int, QString);
        void hideParticles(bool show);
        
        void startTrackingMode();
        void trackingWidgetMoved(QPointF);
        void trackPatch();
        
        void labelOpacityChanged(int value);
        void updateParticles();
        void updateParticles(int n, pi::ParticleVector& particles);

        void showMeanSquares();
        void showCrossCorrelation();
        void showEntropy();
        void showNormalizedEntropy();

        void showGrid(bool);
        
        void run();


        QGraphicsScene* scene(int n);

    protected:
        virtual void mousePressed(EllipseItem* sender, QGraphicsSceneMouseEvent* event);
        virtual void mouseReleased(EllipseItem* sender, QGraphicsSceneMouseEvent* event);

        void mousePressed(QGraphicsRealImageItem*, QGraphicsSceneMouseEvent*);
        void mouseMoved(QGraphicsRealImageItem*, QGraphicsSceneMouseEvent*);
        void mouseReleased(QGraphicsRealImageItem*, QGraphicsSceneMouseEvent*);
        void hoverEntered(QGraphicsRealImageItem*, QGraphicsSceneHoverEvent*);
        void hoverMoved(QGraphicsRealImageItem*, QGraphicsSceneHoverEvent*);
        void hoverLeft(QGraphicsRealImageItem*, QGraphicsSceneHoverEvent*);

    private:
        void connectSignals();
        QGraphicsRealImageItem* showImage(int n, pi::RealImage::Pointer image);
        QGraphicsPixmapItem* showLabel(int n, pi::LabelImage::Pointer image);
        void createParticleItems(int id, int n);
        void rectToRegion(QRectF&, pi::RealImage::RegionType&);

    private:
        friend class ImagePatchInteraction;

        Ui_Simul2D* _ui;
        QWidget* _parent;

        pi::RealImage::Pointer _image[2], _auxImage[2];
        pi::LabelImage::Pointer _label[2];


        pi::PatchTracking _tracking;
        
        QGraphicsScene* _scene[2];
        QGraphicsScene* _miniScene;
        
        QGraphicsRealImageItem* _imageItem[2];
        QGraphicsRealImageItem* _auxImageItem[2];
        QGraphicsRealImageItem* _patchItem[2];

        QGraphicsRectItem* _rects[2];
        QGraphicsRectItem* _selectedRects[2];
        
        QGraphicsPixmapItem* _labelItem[2];
        QVector<EllipseItem*> _particleItem[2];
        QAbstractGraphicsShapeItem* _trackingItem[2];
        QGraphicsRectItem* _trackingRect;
        QGraphicsRectWidget* _trackingWidget;

        pi::ParticleSystemSolver* _solver;

        SelectionMode _particleSelectionMode;
        int _particleSelectedId;

        ImageInteractionMode _imageInteractionMode;
    };
}
#endif /* defined(__ParticleGuidedRegistration__piqSimulCore__) */
