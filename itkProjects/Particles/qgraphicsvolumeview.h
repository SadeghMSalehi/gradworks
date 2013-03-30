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
    enum PixmapItemKey { SliceIndex, RealSliceIndex, AnnotationType };
    enum { SliceImage, SliceMarker, WorkingSet };
    
    QGraphicsVolumeView(QWidget* parent = NULL);

    void setDisplayCollection(pi::AIRDisplayCollection* images);
    void keyReleaseEvent(QKeyEvent* key);
    void clear();

    void mousePressEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    void mouseDoubleClickEvent(QMouseEvent* event);
    
    std::vector<int> getWorkingSet();

signals:
    void sliceDoubleClicked(int n);

public slots:
    void updateDisplay();
    void createWorkingSet();
    void removeWorkingSetItem(int i);
    void clearWorkingSet();
    void currentSliceChanged(int slice);

private:
    bool checkSliceCache();
    bool checkVolumeCache();
    bool updateSource();
    void addWorkingSetItem(QGraphicsPixmapItem *sliceItem);

    
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

    itk::ModifiedTimeType _volumeSourceModifiedTime;
    pi::AIRImage::Pointer _volumeSource;

    typedef QSet<int> QIntSet;
    QIntSet _workingSet;

    QGraphicsRectItem* _currentSliceMarker;
    QVector<QGraphicsPixmapItem*> _slicePixmaps;
};

#endif /* defined(__ParticleGuidedRegistration__qgraphicsvolumeview__) */
