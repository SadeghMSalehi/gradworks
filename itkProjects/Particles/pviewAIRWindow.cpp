//
//  pviewAIRWindow.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/25/13.
//
//

#include "pviewAIRWindow.h"
#include "QFileSystemModel"
#include "QGLWidget"
#include "QVTKWidget2.h"
#include "vtkGenericOpenGLRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkInteractorStyle.h"
#include "vtkBoxWidget.h"
#include "piVTK.h"
#include "vtkProperty.h"
#include "QVTKInteractor.h"
#include "vtkInteractorStyleSwitch.h"
#include "vtkCamera.h"
#include "vnl/vnl_vector.h"
#include "vtkMatrix4x4.h"
#include "piImageIO.h"
#include "itkARGBColorFunction.h"
#include "itkScalarToARGBColormapImageFilter.h"
#include "itkExtractImageFilter.h"
#include "qgraphicscompositeimageitem.h"

using namespace std;
using namespace pi;

typedef itk::Image<RGBAPixel, 3> RGBAVolumeType;

class vtkMouseHandler: public vtkCommand {
private:
    AIRWindow* m_AIRWindow;
    bool m_moving;
public:
    vtkMouseHandler() {
        m_moving = false;
    }

    void SetAIRWindow(AIRWindow* AIRWindow) {
        m_AIRWindow = AIRWindow;
    }

    virtual void Execute(vtkObject *caller, unsigned long eventId, void *callData) {
        vtkInteractorStyle* style = dynamic_cast<vtkInteractorStyle*>(caller);
        if (style == NULL) {
            return;
        }
        vtkRenderWindowInteractor* rwi = style->GetInteractor();
        int* eventPos = rwi->GetEventPosition();
        style->FindPokedRenderer(eventPos[0], eventPos[1]);
        vtkRenderer* currRenderer = style->GetCurrentRenderer();
        if (currRenderer == NULL) {
            return;
        }
        vtkCamera* camera = currRenderer->GetActiveCamera();
        vnl_vector<double> eyePosition(3);
        camera->GetEyePosition(eyePosition.data_block());

        vtkMatrix4x4* viewMat = camera->GetViewTransformMatrix();

        switch (eventId) {
            case vtkCommand::LeftButtonPressEvent: {
                m_moving = true;
                style->OnLeftButtonDown();
                break;
            }
            case vtkCommand::LeftButtonReleaseEvent: {
                m_moving = false;
                style->OnLeftButtonUp();
                break;
            }
            case vtkCommand::MouseMoveEvent: {
                if (m_moving) {
                    m_AIRWindow->UpdateMovingDisplayTransform(viewMat);
                }
                style->OnMouseMove();
                break;
            }
            default:
                return;
        }
    }
};


AIRWindow::AIRWindow(QWidget* parent): m_sliceDirectionActions(this) {
    ui.setupUi(this);
    vtkRenderer* renderer = vtkRenderer::New();
    m_mouseHandler = new vtkMouseHandler();
    m_mouseHandler->SetAIRWindow(this);

    ui.vtkControl->autoBufferSwap();
    ui.vtkControl->autoFillBackground();

    ui.vtkControl->GetRenderWindow()->AddRenderer(renderer);
    vtkInteractorStyle* style = dynamic_cast<vtkInteractorStyle*>(ui.vtkControl->GetInteractor()->GetInteractorStyle());
    if (style != NULL) {
        style->AddObserver(vtkCommand::LeftButtonPressEvent, m_mouseHandler);
        style->AddObserver(vtkCommand::LeftButtonReleaseEvent, m_mouseHandler);
        style->AddObserver(vtkCommand::MouseMoveEvent, m_mouseHandler);
    }
    //    ui.graphicsView->setViewport(qvtkWidget);
    ui.graphicsView->setScene(&m_scene);

    pivtk::PolyDataPointer sphere = pivtk::CreateSphere(60, 60);
    m_propScene.SetRenderer(renderer);
    m_propScene.AddPolyData("sphere", sphere);
    m_propScene.SetRepresentation(VTK_WIREFRAME);

    renderer->ResetCamera();
    renderer->Render();

    m_fixedId = -1;
    m_movingId = -1;

    m_compositeDisplay = new QGraphicsCompositeImageItem();
    m_compositeDisplay->SetImageDisplays(&imageDisplays);
    m_scene.addItem(m_compositeDisplay);

    // Slice Direction
    m_sliceDirectionActions.addAction(ui.actionSliceIJ);
    m_sliceDirectionActions.addAction(ui.actionSliceJK);
    m_sliceDirectionActions.addAction(ui.actionSliceKI);
    m_sliceDirectionActions.setExclusive(true);
    ui.actionSliceIJ->setChecked(true);
    m_currentSliceDir = IJ;
    
    m_compositeDisplay->grabMouse();
    m_compositeDisplay->grabKeyboard();
    m_compositeDisplay->setSelected(true);

    QObject::connect(m_compositeDisplay, SIGNAL(originChanged()), this, SLOT(UpdateOriginDisplay()));
    QObject::connect(ui.actionSliceIJ, SIGNAL(triggered()), this, SLOT(UpdateSliceDirection()));
    QObject::connect(ui.actionSliceJK, SIGNAL(triggered()), this, SLOT(UpdateSliceDirection()));
    QObject::connect(ui.actionSliceKI, SIGNAL(triggered()), this, SLOT(UpdateSliceDirection()));
}

AIRWindow::~AIRWindow() {
    delete m_mouseHandler;
}

void AIRWindow::LoadImage(QString fileName) {
    if (!fileName.isEmpty() && !imageDisplays.IsValidId(m_movingId)) {
        ImageIO<RealImage> io;
        RealImage::Pointer img = io.ReadImage(fileName.toStdString());

        if (imageDisplays.Count() == 0) {
            m_fixedId = 0;
            LoadFixedImage(img);
        } else if (imageDisplays.Count() == 1){
            m_movingId = 1;
            LoadMovingImage(img);
        }
        UpdateCompositeDisplay();

        /*
         // generate some data:
         const int nHistoSize = m_fixedImg.imageHistogram.binData.size();
         QVector<double> x(nHistoSize), y(nHistoSize); // initialize with entries 0..100
         for (int i=0; i<m_fixedImg.imageHistogram.binData.size(); ++i)
         {
         x[i] = i;
         y[i] = m_fixedImg.imageHistogram.binData[i];
         cout << y[i] << endl;
         }
         // create graph and assign data to it:
         ui.histogramView->addGraph();
         ui.histogramView->graph(0)->setData(x, y);
         ui.histogramView->graph(0)->setBrush(QBrush(Qt::yellow));
         // give the axes some labels:
         ui.histogramView->xAxis->setLabel("x");
         ui.histogramView->yAxis->setLabel("y");
         // set axes ranges, so we see all data:
         ui.histogramView->rescaleAxes();
         ui.histogramView->replot();

         ui.intensitySlider->setRealMin(m_fixedImg.imageStat.minIntensity);
         ui.intensitySlider->setRealMax(m_fixedImg.imageStat.maxIntensity);

         ui.intensitySlider->setLowValue(ui.intensitySlider->minimum());
         ui.intensitySlider->setHighValue(ui.intensitySlider->maximum());
         */
    }
}

void AIRWindow::LoadFixedImage(RealImage::Pointer image) {
    ImageDisplayType& dispImg = imageDisplays.AddImage(image);

    imageDisplays.SetReferenceId(m_fixedId);
    
    ui.intensitySlider->setRealMin(dispImg.histogram.dataMin);
    ui.intensitySlider->setRealMax(dispImg.histogram.dataMax);

    ui.intensitySlider->setLowValue(ui.intensitySlider->minimum());
    ui.intensitySlider->setHighValue(ui.intensitySlider->maximum());

    m_currentSliceIndex[IJ] = imageDisplays.GetReferenceSize(IJ)/2.0;
    m_currentSliceIndex[JK] = imageDisplays.GetReferenceSize(JK)/2.0;
    m_currentSliceIndex[KI] = imageDisplays.GetReferenceSize(KI)/2.0;

    // slice axis selection
    UpdateSliceDirection();
}

void AIRWindow::LoadMovingImage(RealImage::Pointer image) {
    ImageDisplayType& dispImg = imageDisplays.AddImage(image);

    ui.ox->setSingleStep(dispImg.srcSpacing[0]);
    ui.oy->setSingleStep(dispImg.srcSpacing[1]);
    ui.oz->setSingleStep(dispImg.srcSpacing[2]);

    ui.intensitySlider2->setRealMin(dispImg.histogram.dataMin);
    ui.intensitySlider2->setRealMax(dispImg.histogram.dataMax);

    ui.intensitySlider2->setLowValue(ui.intensitySlider2->minimum());
    ui.intensitySlider2->setHighValue(ui.intensitySlider2->maximum());

    ui.intensitySlider2->setEnabled(true);

    ui.originBox->setEnabled(true);
    ui.scaleBox->setEnabled(true);

    UpdateOriginDisplay();
}

bool AIRWindow::UpdateMovingDisplayTransform(vtkMatrix4x4* mat) {
    if (m_movingId >= 0 && imageDisplays.Count() > m_movingId) {
        imageDisplays[m_movingId].SetAffineTransform(mat);
        UpdateCompositeDisplay();
        return true;
    }
    return false;
}

void AIRWindow::UpdateCompositeDisplay() {
    m_compositeDisplay->Refresh(m_fixedId, m_movingId);
}

void AIRWindow::on_compositeOpacity_valueChanged(int n) {
    double alpha = 1-double(n)/ui.compositeOpacity->maximum();
    m_compositeDisplay->CompositionModeToAlpha(alpha);
    UpdateCompositeDisplay();
}


void AIRWindow::on_sliceSlider_valueChanged(int n) {
    m_currentSliceIndex[m_currentSliceDir] = n;
    imageDisplays.SetSliceGrid(m_currentSliceDir, m_currentSliceIndex[m_currentSliceDir]);
    UpdateCompositeDisplay();
}

void AIRWindow::on_zoomSlider_valueChanged(int n) {
    QTransform transform;
    transform.scale(n/100.0, n/100.0);
    ui.graphicsView->setTransform(transform);
}

void AIRWindow::on_intensitySlider_lowValueChanged(int n) {
    on_intensitySlider_highValueChanged(0);
}

void AIRWindow::on_intensitySlider_highValueChanged(int n) {
    // correct to handle multiple image histogram intensity
    imageDisplays[m_fixedId].histogram.rangeMin = ui.intensitySlider->realLowValue();
    imageDisplays[m_fixedId].histogram.rangeMax = ui.intensitySlider->realHighValue();
    UpdateCompositeDisplay();
}

void AIRWindow::on_intensitySlider2_lowValueChanged(int n) {
    on_intensitySlider2_highValueChanged(0);
}

void AIRWindow::on_intensitySlider2_highValueChanged(int n) {
    // correct to handle multiple image histogram intensity
    imageDisplays[m_movingId].histogram.rangeMin = ui.intensitySlider2->realLowValue();
    imageDisplays[m_movingId].histogram.rangeMax = ui.intensitySlider2->realHighValue();
    UpdateCompositeDisplay();
}

void AIRWindow::on_alphaOptions_toggled(bool b) {
    if (b) {
        ui.checkerBoardOptions->setChecked(false);
        ui.alphaOptions->setChecked(true);
        ui.compositionOptions->setChecked(false);
        on_compositeOpacity_valueChanged(ui.compositeOpacity->value());
    }
}

void AIRWindow::on_checkerBoardOptions_toggled(bool b) {
    if (b) {
        ui.alphaOptions->setChecked(false);
        ui.checkerBoardOptions->setChecked(true);
        ui.compositionOptions->setChecked(false);
        m_compositeDisplay->CompositionModeToCheckerBoard(ui.checkerBoardRows->value(), ui.checkerBoardCols->value());
        UpdateCompositeDisplay();
    }
}

void AIRWindow::on_compositionOptions_toggled(bool b) {
    if (b) {
        ui.alphaOptions->setChecked(false);
        ui.checkerBoardOptions->setChecked(false);
        ui.compositionOptions->setChecked(true);
        if (ui.compositionDifference->isChecked()) {
            m_compositeDisplay->CompositionMode(0);
        }
        UpdateCompositeDisplay();
    }
}


void AIRWindow::on_ox_valueChanged(double n) {
    if (imageDisplays.IsValidId(m_movingId)) {
        imageDisplays[m_movingId].SetOrigin(0, n);
        UpdateCompositeDisplay();
    }
}

void AIRWindow::on_oy_valueChanged(double n) {
    if (imageDisplays.IsValidId(m_movingId)) {
        imageDisplays[m_movingId].SetOrigin(1, n);
        UpdateCompositeDisplay();
    }
}

void AIRWindow::on_oz_valueChanged(double n) {
    if (imageDisplays.IsValidId(m_movingId)) {
        imageDisplays[m_movingId].SetOrigin(2, n);
        UpdateCompositeDisplay();
    }
}

void AIRWindow::ToggleBlockSignals(bool signalStatus) {
    ui.ox->blockSignals(signalStatus);
    ui.oy->blockSignals(signalStatus);
    ui.oz->blockSignals(signalStatus);
}

void AIRWindow::UpdateOriginDisplay() {
    if (imageDisplays.IsValidId(m_movingId)) {
        ToggleBlockSignals(true);
        ui.ox->setValue(imageDisplays[m_movingId].srcImg->GetOrigin()[0]);
        ui.oy->setValue(imageDisplays[m_movingId].srcImg->GetOrigin()[1]);
        ui.oz->setValue(imageDisplays[m_movingId].srcImg->GetOrigin()[2]);
        ToggleBlockSignals(false);
    }
}

void AIRWindow::UpdateSliceDirection() {
    if (imageDisplays.IsValidId(m_fixedId)) {
        if (ui.actionSliceIJ->isChecked()) {
            m_currentSliceDir = IJ;
        } else if (ui.actionSliceJK->isChecked()) {
            m_currentSliceDir = JK;
        } else if (ui.actionSliceKI->isChecked()) {
            m_currentSliceDir = KI;
        }
        ui.sliceSlider->setValue(m_currentSliceIndex[m_currentSliceDir]);
        ui.sliceSlider->setMaximum(imageDisplays.GetReferenceSize(m_currentSliceDir));
    }
}