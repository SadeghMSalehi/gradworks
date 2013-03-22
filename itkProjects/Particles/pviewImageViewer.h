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


class QGLWidget;
class QVTKWidget2;
class vtkMouseHandler;


class QGraphicsImageItem: public QGraphicsPixmapItem {
public:
    int w;
    int h;

    pi::RealImage::Pointer getImage() { return _srcImg; }
    void setWindowRange(pi::DataReal min, pi::DataReal max) {
        _windowMin = min;
        _windowMax = max;
        updatePixmap();
    }
    void setImage(pi::RealImage::Pointer srcImg);

private:
    void updatePixmap();

private:
    pi::DataReal _windowMin, _windowMax;
    pi::RealImage::Pointer _srcImg;
};

class ImageViewer : public QDialog {
    Q_OBJECT
public:
    typedef itk::Image<pi::DataReal,2> RealImage2D;
    pi::RealImage::Pointer fixedImg;

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

private:
    friend class vtkMouseHandler;
    
    Ui::ImageViewer ui;
    QGraphicsScene m_scene;
    QVTKWidget2* qvtkWidget;
    pivtk::PropScene m_propScene;
    pi::ImageReslice<pi::RealImage> m_fixedImg;
    pi::ImageReslice<pi::RealImage> m_movingImg;
    QGraphicsPixmapItem* m_fixedPixmap;
    QGraphicsPixmapItem* m_movingPixmap;
    vtkMouseHandler* m_mouseHandler;
};

#endif /* defined(__ParticleGuidedRegistration__pviewImageTransform__) */
