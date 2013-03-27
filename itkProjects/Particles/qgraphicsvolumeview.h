//
//  qgraphicsvolumeview.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/27/13.
//
//

#ifndef __ParticleGuidedRegistration__qgraphicsvolumeview__
#define __ParticleGuidedRegistration__qgraphicsvolumeview__

#include <iostream>

#include <QGraphicsView>
#include <QGraphicsScene>

#include "piImageSlice.h"

class QGraphicsVolumeView: public QGraphicsView {
    Q_OBJECT
public:
    QGraphicsVolumeView(QWidget* parent = NULL);

    void setDisplayCollection(pi::AIRDisplayCollection* images);

public slots:
    void updateDisplay();

private:
    bool checkSliceCache();
    bool checkVolumeCache();
    bool updateSource();

private:
    int _thumbsWidth;
    int _columnCount;

    int _displayId;
    bool _displayReference;
    bool _manualIntensityScaling;

    QGraphicsScene _scene;
    pi::AIRDisplayCollection* _airImages;
    pi::AIRImageVector _sliceCache;
    pi::RGBAImageVector _displayImages;
    pi::SliceDirectionEnum _directionCache;
    pi::AIRImage::Pointer _volumeCache;
};

#endif /* defined(__ParticleGuidedRegistration__qgraphicsvolumeview__) */
