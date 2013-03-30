//
//  qgraphicslablemapitem.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/29/13.
//
//

#ifndef __ParticleGuidedRegistration__qgraphicslablemapitem__
#define __ParticleGuidedRegistration__qgraphicslablemapitem__

#include <iostream>
#include <QGraphicsPixmapItem>
#include "piImageSlice.h"

class QGraphicsLabelMapItem: public QObject, public QGraphicsPixmapItem {
    Q_OBJECT
private:
    pi::AIRLabelSlice _volumeSlice;
    
public:
    QGraphicsLabelMapItem(QGraphicsItem* parent);
    ~QGraphicsLabelMapItem();
    
    void setLabelVolume();
    void setColorMap(uint* colorMap);
    void changeDirection(pi::SliceDirectionEnum dir);
    void changeSlice(int slice);
    
    pi::AIRLabelSlice& getSlice();

private:
    void sliceToVolume();
    
public slots:
    void labelVolumeRequested();
    
};

#endif /* defined(__ParticleGuidedRegistration__qgraphicslablemapitem__) */
