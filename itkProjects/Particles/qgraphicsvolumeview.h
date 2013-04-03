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
    virtual ~QGraphicsVolumeView() {
        if (_scene != NULL) {
            delete _scene;
        }
    }

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
    void directionChanged(pi::SliceDirectionEnum dir);

    void moveToVolume(int i);
    void flipLR(bool toggle);
    void flipUD(bool toggle);

private:
    void addWorkingSetItem(QGraphicsItem *sliceItem);
    void updatePixmaps();
    
private:
    // volume view storage
    typedef pi::VolumeDisplay<pi::AIRImage> AIRVolumeDisplay;
    QVector<AIRVolumeDisplay> _volumeDisplays;

//    typedef QList<QGraphicsItem*> QGraphicsItemList;
//    typedef QList<QGraphicsItemList> QGraphicsItemTable;
//    QGraphicsItemTable _pixmapTable;
//    QGraphicsItemTable _pximapRegionTable;

    bool _xFlipped;
    bool _yFlipped;
    bool _showAll;

    int _thumbsWidth;
    double _rescaleFactor;
    int _columnCount;
    bool _useNavigationImage;
    pi::SliceDirectionEnum _directionCache;

    QGraphicsScene* _scene;
    pi::AIRDisplayCollection* _airImages;

    typedef QSet<int> QIntSet;
    QIntSet _workingSet;

    QGraphicsRectItem* _currentSliceMarker;
};

#endif /* defined(__ParticleGuidedRegistration__qgraphicsvolumeview__) */
