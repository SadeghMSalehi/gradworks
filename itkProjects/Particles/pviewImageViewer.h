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


class QGraphicsVolumeItem: public QGraphicsPixmapItem {
public:
    typedef itk::AffineTransform<double> TransformType;

    int w;
    int h;

    QGraphicsVolumeItem(QGraphicsItem* parent = NULL);
    
    pi::RealImage::Pointer getImage() { return _srcImg; }
    void setWindowRange(pi::DataReal min, pi::DataReal max) {
        _windowMin = min;
        _windowMax = max;
        updatePixmap();
    }
    void setImage(pi::RealImage::Pointer srcImg);
    void setResampleGrid(pi::RealImage::Pointer grid);
    void setTransform(vtkMatrix4x4* mat);

private:
    void resampleGrid();
    void updatePixmap();

private:
    pi::DataReal _windowMin, _windowMax;
    pi::RealImage::Pointer _srcImg, _resampleGrid, _resampledImg;
    TransformType::Pointer _transform;
};


class ImageViewer : public QDialog {
    Q_OBJECT
public:
    typedef itk::Image<pi::DataReal,2> RealImage2D;
    pi::RealImage::Pointer fixedImg, fixedGrid, movingImg, movingGrid;

public:
    ImageViewer(QWidget* parent = NULL);
    virtual ~ImageViewer();
    void LoadImage(QString fileName);
    void LoadMovingImage(QString fileName);
    void ResampleSlice();

public slots:
    void on_fixedOpacity_sliderMoved(int n);
    void on_movingOpacity_sliderMoved(int n);

    void on_fixedSliceSlider_sliderMoved(int n);
    void on_zoomSlider_sliderMoved(int n);

    void on_intensitySlider_lowValueChanged(int n);
    void on_intensitySlider_highValueChanged(int n);
    void on_intensitySlider2_lowValueChanged(int n);
    void on_intensitySlider2_highValueChanged(int n);

private:
    void SetFixedSlice(int dir, int idx = -1);
    void SetMovingSlice(int dir, int idx = -1);

private:
    friend class vtkMouseHandler;
    
    Ui::ImageViewer ui;
    QGraphicsScene m_scene;
    QVTKWidget2* qvtkWidget;
    pivtk::PropScene m_propScene;
    QGraphicsVolumeItem* m_fixedItem;
    QGraphicsVolumeItem* m_movingItem;
    vtkMouseHandler* m_mouseHandler;
};

#endif /* defined(__ParticleGuidedRegistration__pviewImageTransform__) */
