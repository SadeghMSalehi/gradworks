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
#include "piVolumeDisplay.h"


class QGraphicsVolumeView: public QGraphicsView {
    Q_OBJECT
public:
    enum PixmapItemKey { SliceIndex, RealSliceIndex, AnnotationType };
    enum { SliceImage, SliceMarker, WorkingSet };
    
    QGraphicsVolumeView(QWidget* parent = NULL);
    virtual ~QGraphicsVolumeView() {}

    void setThumbsWidth(int w) { _thumbsWidth = w; }
    void setDisplayCollection(pi::AIRDisplayCollection* images, bool useNavigationImage = true);
    void keyReleaseEvent(QKeyEvent* key);
    void clear();
    void fitToImage(int sliceIdx, int volume = -1);
    std::vector<int> getWorkingSet();

protected:
    virtual void mousePressEvent(QMouseEvent* event);
    virtual void mouseReleaseEvent(QMouseEvent* event);
    virtual void mouseDoubleClickEvent(QMouseEvent* event);

    virtual void dragEnterEvent(QDragEnterEvent *event);
    virtual void dropEvent(QDropEvent *event);

signals:
    void sliceDoubleClicked(int n);
    void fileDropped(QString fileName);

public slots:
    void updateDisplay(int volumeId = -1);
    void createWorkingSet();
    void removeWorkingSetItem(int i);
    void clearWorkingSet();
    void currentSliceChanged(int slice);
    void setVolumeToShow(int i);
    void moveToVolume(int i);
    void setImageFlip(bool xFlip, bool yFlip);

private:
    void addWorkingSetItem(QGraphicsPixmapItem *sliceItem);

    
private:
    // volume view storage
    typedef pi::VolumeDisplay<pi::AIRImage> AIRVolumeDisplay;
    typedef QHash<int, AIRVolumeDisplay> VolumeHash;
    VolumeHash _volumeDisplays;

    bool _xFlipped;
    bool _yFlipped;

    int _thumbsWidth;
    double _rescaleFactor;
    int _columnCount;
//
//    int _displayId;
//    bool _displayReference;
//    bool _manualIntensityScaling;
    bool _useNavigationImage;
    pi::SliceDirectionEnum _directionCache;


    QGraphicsScene _scene;
    pi::AIRDisplayCollection* _airImages;
//    pi::AIRImageVector _sliceCache;
//    pi::RGBAImageVector _displayImages;
//    pi::AIRImage::Pointer _volumeCache;
//    itk::ModifiedTimeType _volumeSourceModifiedTime;
//    pi::AIRImage::Pointer _volumeSource;

    typedef QSet<int> QIntSet;
    QIntSet _workingSet;

    QGraphicsRectItem* _currentSliceMarker;
//    QVector<QGraphicsPixmapItem*> _slicePixmaps;
};

#endif /* defined(__ParticleGuidedRegistration__qgraphicsvolumeview__) */
