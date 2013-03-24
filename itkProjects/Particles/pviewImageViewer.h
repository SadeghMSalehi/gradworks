//
//  pviewImageTransform.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/14/13.
//
//

#ifndef __ParticleGuidedRegistration__pviewImageTransform__
#define __ParticleGuidedRegistration__pviewImageTransform__

#include <iostream>
#include "ui_pviewImageViewer.h"
#include "QGraphicsScene"
#include "QGraphicsPixmapItem"
#include "piVTK.h"
#include "piImageDef.h"
#include "piImageSlice.h"
#include "vtkMatrix4x4.h"


class QGLWidget;
class QVTKWidget2;
class vtkMouseHandler;

typedef typename pi::ImageDisplayCollection<pi::RealImage>::ImageDisplayType ImageDisplayType;

class QGraphicsCompositeImageItem: public QGraphicsPixmapItem {
private:
    enum CompositionMode { Alpha, CheckerBoard };

private:
    int _resampleIdx;
    double _alpha;
    int _cbRows, _cbCols;
    CompositionMode _compositionMode;


    // memory holder for qpixmap
    pi::RGBAVolumeType::Pointer _rgbImage;

    // memory holder for composite image
    pi::RealImage::Pointer _compositeImage;

private:
    typename pi::ImageDisplayCollection<pi::RealImage>* _imageDisplays;

    int _fixedId;
    int _movingId;

    // image displays

    void CompositeAlpha();
    
    inline pi::DataReal Clamp(pi::DataReal f, const pi::DataReal fMin, const pi::DataReal fMax) {
        if (f < fMin) {
            f = fMin;
        } else if (f > fMax) {
            f = fMax;
        }
        f = (f-fMin)/(fMax-fMin)*65535;
        return f;
    }

public:
    QGraphicsCompositeImageItem(QGraphicsItem* parent = NULL): QGraphicsPixmapItem(parent), _resampleIdx(0) {
        _fixedId = 0;
        _movingId = 1;
        _imageDisplays = NULL;
        _compositionMode = Alpha;
        _alpha = 1;
        _cbRows = 4;
        _cbCols = 4;
    }

    void SetImageDisplays(typename pi::ImageDisplayCollection<pi::RealImage>* displays) {
        _imageDisplays = displays;
    }

    void CompositionModeToAlpha(double alpha) {
        _compositionMode = Alpha;
        _alpha = alpha;
    }

    void CompositionModeToCheckerBoard(int r, int c) {
        _compositionMode = CheckerBoard;
        _cbRows = r;
        _cbCols = c;
    }

    void Refresh();

};

//
//class QGraphicsVolumeItem: public QGraphicsPixmapItem {
//public:
//    enum CompositionMode { DIFF, XOR };
//    pi::ImageDisplayCollection<pi::RealImage> _imageDisplays;
//
//    QGraphicsVolumeItem(QGraphicsItem* parent = NULL);
//
//    void setResampleGrid(pi::RealImage::Pointer grid) {
//        _resampleGrid = grid;
//    }
//
//    void setTransform(int n, vtkMatrix4x4* mat) {
//        _imageDisplays[n].SetAffineTransform(mat);
//    }
//
//    void setAlphaBlending(int o1, int o2, double alpha);
//    void setComposition(int o1, int o2, CompositionMode mode);
//    void setCheckerBoard(int o1, int o2, int r, int c);
//    void setDisplay(int o);
//
//private:
//    void updatePixmap();
//
//private:
//    pi::RealImage::Pointer _resampleGrid;
//    pi::RealImage::Pointer _displayImage;
//    pi::RGBAVolumeType::Pointer _rgbDisplay;
//    std::vector<ImageDisplayType> _imageDisplays;
//
//    int _pixmapWidth;
//    int _pixmapHeight;
//};


class ImageViewer : public QDialog {
    Q_OBJECT
public:
    pi::ImageDisplayCollection<pi::RealImage> imageDisplays;
    
public:
    ImageViewer(QWidget* parent = NULL);
    virtual ~ImageViewer();
    
    void LoadImage(QString fileName);
    void ClearImages();

public slots:
    void on_compositeOpacity_sliderMoved(int n);

    void on_fixedSliceSlider_sliderMoved(int n);
    void on_zoomSlider_sliderMoved(int n);

    void on_intensitySlider_lowValueChanged(int n);
    void on_intensitySlider_highValueChanged(int n);


private:
    friend class vtkMouseHandler;
    
    Ui::ImageViewer ui;
    QGraphicsScene m_scene;
    QVTKWidget2* qvtkWidget;
    pivtk::PropScene m_propScene;
    vtkMouseHandler* m_mouseHandler;


    QGraphicsCompositeImageItem* m_compositeDisplay;
};

#endif /* defined(__ParticleGuidedRegistration__pviewImageTransform__) */
