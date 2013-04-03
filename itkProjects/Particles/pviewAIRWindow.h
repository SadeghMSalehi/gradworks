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
#include "piImageIO.h"
#include "vtkMatrix4x4.h"

namespace air {
    class AlgorithmManager;
}

class QGLWidget;
class QVTKWidget2;
class vtkMouseHandler;
class QGraphicsCompositeImageItem;

typedef typename pi::ImageDisplayCollection<pi::AIRImage>::ImageDisplayType ImageDisplayType;
typedef typename pi::ImageDisplayCollection<pi::AIRImage>::ImageDisplayPointer ImageDisplayPointer;

class AIRWindow : public QMainWindow {
    Q_OBJECT
public:
    pi::ImageDisplayCollection<pi::AIRImage> imageDisplays;

signals:
    void fileDropped(QString& fileName);

public:
    AIRWindow(QWidget* parent = NULL);
    virtual ~AIRWindow();

    void LoadImage(QString fileName);
    void ClearImages();

    void UpdateCompositeDisplay();
    bool UpdateMovingDisplayTransform(vtkMatrix4x4* mat);

    bool IsImage1Loaded();
    bool IsImage2Loaded();

    int getForegroundLabel();
    int getBackgroundLabel();
    
protected:
    void dragEnterEvent(QDragEnterEvent* event);
    void dropEvent(QDropEvent* event);

    void setupShortcutKeys();
    
public slots:
    void on_actionTest_triggered();
    void on_actionResample_triggered();
    void on_actionLoadTransform_triggered();
    void on_actionSaveTransform_triggered();
    void on_actionMultipleSlice_triggered();
    void on_actionUnload_triggered();

    void on_compositeOpacity_valueChanged(int n);

    void on_sliceSlider_valueChanged(int n);
    void on_zoomSlider_valueChanged(int n);

    void on_intensitySlider_lowValueChanged(int n);
    void on_intensitySlider_highValueChanged(int n);
    void on_intensitySlider_sliderMoved(int n);
    void on_intensitySlider2_lowValueChanged(int n);
    void on_intensitySlider2_highValueChanged(int n);

    void on_checkerBoardOptions_toggled(bool checked);
    void on_alphaOptions_toggled(bool checked);
    void on_compositionOptions_toggled(bool checked);

    void OnTranslationWidgetChanged();
    void UpdateTranslationWidget();
    void ChangeSliceDirection();
    void ChangeSliceIndex(int n);
    void ChangeInteractionMode();
    void LoadSegmentation();
    void SaveSegmentation();
    void SaveSegmentationAs();
    void PropagateSegmentation();
    void brushLabelChanged(int n);
    
    void on_image1Name_clicked(bool checked);
    void on_image2Name_clicked(bool checked);
    void on_image1Name_fileDropped(QString& fileName);
    void on_image2Name_fileDropped(QString& fileName);

    void previousSlice();
    void nextSlice();
    void sliceZoomIn();
    void sliceZoomOut();
    void labelOpacityUp();
    void labelOpacityDown();
    void copyLabel();


    // test
    void showTestDialog();
    
private:
    friend class vtkMouseHandler;
    friend class air::AlgorithmManager;
    
    Ui::AIRWindow ui;
    QGraphicsScene m_scene;
    QVTKWidget2* qvtkWidget;
    pivtk::PropScene m_propScene;
    vtkMouseHandler* m_mouseHandler;

    int m_fixedId;
    int m_movingId;
    pi::SliceDirectionEnum m_requestedSliceDir, m_currentSliceDir;
    int m_currentSliceIndex[3];

    QActionGroup m_sliceDirectionActions;
    QGraphicsCompositeImageItem* m_compositeDisplay;

    pi::ImageIO<pi::AIRImage> io;
    QButtonGroup* _drawingToolButtons;

    QList<QIcon> _segmentationLabels;

    air::AlgorithmManager* _algorithms;
};

#endif /* defined(__ParticleGuidedRegistration__pviewAIRWindow__) */
