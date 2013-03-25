//
//  pviewAIRWindow.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/25/13.
//
//

#ifndef __ParticleGuidedRegistration__pviewAIRWindow__
#define __ParticleGuidedRegistration__pviewAIRWindow__

#include <iostream>
#include "ui_pviewAIRWindow.h"
#include "QGraphicsScene"
#include "QGraphicsPixmapItem"
#include "QActionGroup"

#include "piVTK.h"
#include "piImageDef.h"
#include "piImageSlice.h"
#include "vtkMatrix4x4.h"


class QGLWidget;
class QVTKWidget2;
class vtkMouseHandler;
class QGraphicsCompositeImageItem;

typedef typename pi::ImageDisplayCollection<pi::RealImage>::ImageDisplayType ImageDisplayType;

class AIRWindow : public QMainWindow {
    Q_OBJECT
public:
    pi::ImageDisplayCollection<pi::RealImage> imageDisplays;

public:
    AIRWindow(QWidget* parent = NULL);
    virtual ~AIRWindow();

    void LoadImage(QString fileName);
    void ClearImages();

    void UpdateCompositeDisplay();
    bool UpdateMovingDisplayTransform(vtkMatrix4x4* mat);

public slots:
    void on_compositeOpacity_valueChanged(int n);

    void on_sliceSlider_valueChanged(int n);
    void on_zoomSlider_valueChanged(int n);

    void on_intensitySlider_lowValueChanged(int n);
    void on_intensitySlider_highValueChanged(int n);
    void on_intensitySlider2_lowValueChanged(int n);
    void on_intensitySlider2_highValueChanged(int n);

    void on_checkerBoardOptions_toggled(bool checked);
    void on_alphaOptions_toggled(bool checked);
    void on_compositionOptions_toggled(bool checked);

    void on_ox_valueChanged(double n);
    void on_oy_valueChanged(double n);
    void on_oz_valueChanged(double n);

    void UpdateOriginDisplay();
    void UpdateSliceDirection();

private:
    void LoadFixedImage(pi::RealImage::Pointer image);
    void LoadMovingImage(pi::RealImage::Pointer image);
    void ToggleBlockSignals(bool signalStatus);
private:
    friend class vtkMouseHandler;

    Ui::AIRWindow ui;
    QGraphicsScene m_scene;
    QVTKWidget2* qvtkWidget;
    pivtk::PropScene m_propScene;
    vtkMouseHandler* m_mouseHandler;

    int m_fixedId;
    int m_movingId;
    pi::SliceDirectionEnum m_currentSliceDir;
    int m_currentSliceIndex[3];

    QActionGroup m_sliceDirectionActions;
    QGraphicsCompositeImageItem* m_compositeDisplay;
};

#endif /* defined(__ParticleGuidedRegistration__pviewAIRWindow__) */
