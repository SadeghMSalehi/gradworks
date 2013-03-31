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

    void setThumbsWidth(int w) { _thumbsWidth = w; }
    void setDisplayCollection(pi::AIRDisplayCollection* images, bool useNavigationImage = true);
    void keyReleaseEvent(QKeyEvent* key);
    void clear();
    void fitToImage(int sliceIdx, int volume = 0);

protected:
    void mousePressEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    void mouseDoubleClickEvent(QMouseEvent* event);

    void dragEnterEvent(QDragEnterEvent *event);
    void dropEvent(QDropEvent *event);

    std::vector<int> getWorkingSet();

signals:
    void sliceDoubleClicked(int n);
    void fileDropped(QString fileName);

public slots:
    void updateDisplay();
    void createWorkingSet();
    void removeWorkingSetItem(int i);
    void clearWorkingSet();
    void currentSliceChanged(int slice);
    void setVolumeToShow(int i);

private:
    bool checkSliceCache();
    bool checkVolumeCache();
    bool updateSource();
    void addWorkingSetItem(QGraphicsPixmapItem *sliceItem);

    
private:
    // volume view storage
    typedef pi::VolumeDisplay<pi::AIRImage> AIRVolumeDisplay;
    typedef QHash<int, AIRVolumeDisplay> VolumeHash;
    VolumeHash _volumeDisplays;

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
