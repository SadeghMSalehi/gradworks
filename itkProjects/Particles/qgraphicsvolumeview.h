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
#include <vector>

#include <QGraphicsView>
#include <QGraphicsScene>
#include <QSet>
#include <QVector>

#include "piImageSlice.h"

class QGraphicsVolumeView: public QGraphicsView {
    Q_OBJECT
public:
    enum PixmapItemKey { SliceIndex, SliceHighlight };
    
    QGraphicsVolumeView(QWidget* parent = NULL);

    void setDisplayCollection(pi::AIRDisplayCollection* images);
    void keyReleaseEvent(QKeyEvent* key);
    void clear();

    void mousePressEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    
    std::vector<int> getWorkingSet();

public slots:
    void updateDisplay();
    void createWorkingSet();
    void clearWorkingSet();

private:
    bool checkSliceCache();
    bool checkVolumeCache();
    bool updateSource();
    void highlightSlice(QGraphicsPixmapItem *sliceItem);

    
private:
    int _thumbsWidth;
    double _rescaleFactor;
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

    typedef QSet<int> QIntSet;
    QIntSet _workingSet;

    QVector<QGraphicsPixmapItem*> _slicePixmaps;
};

#endif /* defined(__ParticleGuidedRegistration__qgraphicsvolumeview__) */
