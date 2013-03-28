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
#include "itkARGBColorFunction.h"
#include "itkScalarToARGBColormapImageFilter.h"
#include "itkExtractImageFilter.h"
#include "qgraphicscompositeimageitem.h"
#include "qgraphicsvolumeview.h"
#include "QFileDialog"
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QMimeData>
#include <QList>
#include <QUrl>

typedef QList<QUrl> UrlList;

using namespace std;
using namespace pi;

typedef itk::Image<RGBAPixel, 3> RGBAVolumeType;

class vtkMouseHandler: public vtkCommand {
private:
    AIRWindow* m_AIRWindow;
    bool m_moving;
public:
    vtkTypeMacro(vtkMouseHandler, vtkCommand);
    
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

    // Slice Direction
    m_sliceDirectionActions.addAction(ui.actionSliceIJ);
    m_sliceDirectionActions.addAction(ui.actionSliceJK);
    m_sliceDirectionActions.addAction(ui.actionSliceKI);
    m_sliceDirectionActions.setExclusive(true);
    ui.actionSliceIJ->setChecked(true);
    m_currentSliceDir = IJ;


    connect(this, SIGNAL(fileDropped(QString&)), this, SLOT(on_image1Name_fileDropped(QString&)));

    connect(m_compositeDisplay, SIGNAL(translationChanged()), this, SLOT(UpdateTranslationWidget()));
    connect(ui.actionSliceIJ, SIGNAL(triggered()), this, SLOT(UpdateSliceDirection()));
    connect(ui.actionSliceJK, SIGNAL(triggered()), this, SLOT(UpdateSliceDirection()));
    connect(ui.actionSliceKI, SIGNAL(triggered()), this, SLOT(UpdateSliceDirection()));
    connect(ui.actionDrawing, SIGNAL(triggered()), this, SLOT(ChangeInteractionMode()));

    connect(ui.actionNewWorkingSet, SIGNAL(triggered()), ui.multipleSliceView, SLOT(createWorkingSet()));
    connect(ui.actionClearWorkingSet, SIGNAL(triggered()), ui.multipleSliceView, SLOT(clearWorkingSet()));

    ui.multipleSliceView->addAction(ui.actionNewWorkingSet);
    ui.multipleSliceView->addAction(ui.actionPropagateLabel);
    ui.multipleSliceView->addAction(ui.actionClearWorkingSet);

    ui.multipleSliceView->hide();
}

AIRWindow::~AIRWindow() {
    m_mouseHandler->UnRegister();
}

void AIRWindow::LoadImage(QString fileName) {
    if (!fileName.isEmpty() && !imageDisplays.IsValidId(m_movingId)) {

        if (imageDisplays.Count() == 0) {
            on_image1Name_fileDropped(fileName);
        } else if (imageDisplays.Count() == 1){
            on_image2Name_fileDropped(fileName);
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


void AIRWindow::UpdateCompositeDisplay() {
    m_compositeDisplay->Refresh(m_fixedId, m_movingId);
}


bool AIRWindow::UpdateMovingDisplayTransform(vtkMatrix4x4* mat) {
    if (m_movingId >= 0 && imageDisplays.Count() > m_movingId) {
        imageDisplays[m_movingId].SetAffineTransform(mat);
        UpdateCompositeDisplay();
        return true;
    }
    return false;
}



void AIRWindow::dragEnterEvent(QDragEnterEvent *event) {
    if (event->mimeData()->hasUrls()) {
        UrlList urls = event->mimeData()->urls();
        UrlList::ConstIterator iter;
        for (iter = urls.constBegin(); iter != urls.constEnd(); iter++) {
            const QUrl& url = (*iter);
            if (url.scheme() != "file") {
                event->setAccepted(false);
                return;
            }
        }
        event->setAccepted(true);
        event->acceptProposedAction();
    }
}

void AIRWindow::dropEvent(QDropEvent *event) {
    QList<QString> fileNames;
    if (event->mimeData()->hasUrls()) {
        UrlList urls = event->mimeData()->urls();
        const QUrl& url = urls[0];
        if (url.scheme() != "file") {
            return;
        }
        event->acceptProposedAction();
        QString filePath = url.path();
        emit fileDropped(filePath);
    }
}

#pragma mark -
#pragma mark SINGAL HANDLERS

void AIRWindow::on_image1Name_clicked(bool checked) {
    if (checked) {
        QString fileName = QFileDialog::getOpenFileName(this, "Load a fixed image", ".");
        if (fileName.isNull()) {
            return;
        }
        on_image1Name_fileDropped(fileName);
        if (!imageDisplays.IsValidId(m_fixedId)) {
            ui.image1Name->setChecked(false);
        }
    }
}

void AIRWindow::on_image2Name_clicked(bool checked) {
    if (checked) {
        QString fileName = QFileDialog::getOpenFileName(this, "Load a fixed image", ".");
        if (fileName.isNull()) {
            return;
        }
        on_image2Name_fileDropped(fileName);
        if (!imageDisplays.IsValidId(m_movingId)) {
            ui.image2Name->setChecked(false);
        }
    }
}


void AIRWindow::on_image1Name_fileDropped(QString& fileName) {
    AIRImage::Pointer image = io.ReadCastedImage(fileName.toStdString());
    if (image.IsNull()) {
        return;
    }

    if (!imageDisplays.SetImage(0, image)) {
        return;
    }

    m_fixedId = 0;
    ui.alphaOptions->setEnabled(true);
    ui.compositionOptions->setEnabled(false);
    ui.checkerBoardOptions->setEnabled(false);
    ui.sliceSlider->setEnabled(true);
    ui.intensitySlider->setEnabled(true);
    m_sliceDirectionActions.setEnabled(true);
    ui.actionSliceIJ->setChecked(true);
    ui.actionDrawing->setEnabled(true);
    ui.actionDrawing->setChecked(false);


    ImageDisplayType& dispImg = imageDisplays.GetLast();
    imageDisplays.SetReferenceId(m_fixedId);

    m_scene.addItem(m_compositeDisplay);
    m_compositeDisplay->grabMouse();
    m_compositeDisplay->grabKeyboard();
    m_compositeDisplay->setSelected(true);

    ui.intensitySlider->setRealMin(dispImg.histogram.dataMin);
    ui.intensitySlider->setRealMax(dispImg.histogram.dataMax);

    ui.intensitySlider->setLowValue(ui.intensitySlider->minimum());
    ui.intensitySlider->setHighValue(ui.intensitySlider->maximum());

    m_currentSliceIndex[IJ] = imageDisplays.GetReferenceSize(IJ)/2.0;
    m_currentSliceIndex[JK] = imageDisplays.GetReferenceSize(JK)/2.0;
    m_currentSliceIndex[KI] = imageDisplays.GetReferenceSize(KI)/2.0;

    // slice axis selection
    UpdateSliceDirection();

    ui.image1Name->setText(fileName.right(20));
    ui.image1Name->setToolTip(fileName);
    ui.image1Name->setChecked(true);

    ui.graphicsView->fitInView(m_compositeDisplay->boundingRect(), Qt::KeepAspectRatio);
}

void AIRWindow::on_image2Name_fileDropped(QString& fileName) {
    AIRImage::Pointer image = io.ReadCastedImage(fileName.toStdString());
    if (image.IsNull()) {
        return;
    }

    if (!imageDisplays.SetImage(1, image)) {
        return;
    }

    m_movingId = 1;
    ImageDisplayType& dispImg = imageDisplays.GetLast();
    ImageDisplayType::VectorType tx = dispImg.GetAffineTranslation();
    ui.ox->setSingleStep(tx[0]);
    ui.oy->setSingleStep(tx[1]);
    ui.oz->setSingleStep(tx[2]);

    ui.intensitySlider2->setRealMin(dispImg.histogram.dataMin);
    ui.intensitySlider2->setRealMax(dispImg.histogram.dataMax);

    ui.intensitySlider2->setLowValue(ui.intensitySlider2->minimum());
    ui.intensitySlider2->setHighValue(ui.intensitySlider2->maximum());

    ui.intensitySlider2->setEnabled(true);

    ui.originBox->setEnabled(true);
    ui.scaleBox->setEnabled(true);

    ui.compositionOptions->setEnabled(true);
    ui.checkerBoardOptions->setEnabled(true);

    OnTranslationWidgetChanged();

    ui.image2Name->setText(fileName.right(20));
    ui.image2Name->setToolTip(fileName);
    ui.image2Name->setChecked(true);

}

void AIRWindow::on_actionDrawing_triggered(bool drawing) {
    
}

void AIRWindow::on_actionResample_triggered() {
    if (imageDisplays.IsValidId(m_movingId)) {
        QString saveFilename = QFileDialog::getSaveFileName(this, "Save Image File", ".");
        if (saveFilename.isNull()) {
            return;
        }
        pi::ImageIO<AIRImage> io;
        AIRImage::Pointer resampledImg = imageDisplays[m_movingId].Resample3D(imageDisplays.GetReferenceGrid(), 0);
        io.WriteImage(saveFilename.toUtf8().data(), resampledImg);
    }
}

void AIRWindow::on_actionLoadTransform_triggered() {
    if (imageDisplays.IsValidId(m_movingId)) {
        QString loadFilename = QFileDialog::getOpenFileName(this, "Load Transform File", ".");
        if (loadFilename.isNull()) {
            return;
        }
        pi::ImageIO<AIRImage> io;
        pi::ImageIO<AIRImage>::TransformType::Pointer transform = io.ReadTransform(loadFilename.toUtf8().data());
        imageDisplays[m_movingId].SetAffineTransform(transform);
        UpdateTranslationWidget();
        UpdateCompositeDisplay();
    }
}

void AIRWindow::on_actionSaveTransform_triggered() {
    if (imageDisplays.IsValidId(m_movingId)) {
        QString saveFilename = QFileDialog::getSaveFileName(this, "Save Transform File", ".");
        if (saveFilename.isNull()) {
            return;
        }

        itk::TransformFileWriter::Pointer writer = itk::TransformFileWriter::New();
        writer->SetInput(imageDisplays[m_movingId].GetAffineTransform());
        writer->SetFileName(saveFilename.toUtf8().data());
        writer->Update();
    }
}

void AIRWindow::on_actionMultipleSlice_triggered() {
    if (ui.multipleSliceView->isHidden()) {
        ui.multipleSliceView->show();
    } else {
        ui.multipleSliceView->hide();
    }
    if (imageDisplays.IsValidId(m_fixedId)) {
        if (!ui.multipleSliceView->isHidden()) {
            ui.multipleSliceView->setDisplayCollection(&imageDisplays);
            ui.multipleSliceView->updateDisplay();
        }
    }
}

void AIRWindow::on_actionUnload_triggered() {
    imageDisplays.Reset();
    m_fixedId = -1;
    m_movingId = -1;
    ui.alphaOptions->setEnabled(false);
    ui.compositionOptions->setEnabled(false);
    ui.checkerBoardOptions->setEnabled(false);
    ui.sliceSlider->setEnabled(false);
    ui.sliceSlider->setValue(0);
    ui.intensitySlider->setEnabled(false);
    ui.intensitySlider2->setEnabled(false);
    m_sliceDirectionActions.setEnabled(false);
    ui.actionDrawing->setEnabled(false);
    ui.image1Name->setChecked(false);
    ui.image1Name->setText("Not Loaded");
    ui.image2Name->setText("Not Loaded");
    ui.image1Name->setToolTip("");
    ui.image2Name->setToolTip("");
    ui.image2Name->setChecked(false);
    m_scene.removeItem(m_compositeDisplay);
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

void AIRWindow::on_intensitySlider_sliderMoved(int n) {
    ui.multipleSliceView->updateDisplay();
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


void AIRWindow::OnTranslationWidgetChanged() {
    if (imageDisplays.IsValidId(m_movingId)) {
        ImageDisplayType::VectorType translation;
        translation[0] = ui.ox->value();
        translation[1] = ui.oy->value();
        translation[2] = ui.oz->value();
        imageDisplays[m_movingId].SetAffineTranslation(translation);
        UpdateCompositeDisplay();
    }
}

void AIRWindow::UpdateTranslationWidget() {
    if (imageDisplays.IsValidId(m_movingId)) {
        ui.ox->blockSignals(true);
        ui.oy->blockSignals(true);
        ui.oz->blockSignals(true);
        ImageDisplayType::VectorType translation = imageDisplays[m_movingId].GetAffineTranslation();
        ui.ox->setValue(translation[0]);
        ui.oy->setValue(translation[1]);
        ui.oz->setValue(translation[2]);
        ui.ox->blockSignals(false);
        ui.oy->blockSignals(false);
        ui.oz->blockSignals(false);
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
        imageDisplays.SetSliceGrid(m_currentSliceDir, m_currentSliceIndex[m_currentSliceDir]);
        ui.sliceSlider->setValue(m_currentSliceIndex[m_currentSliceDir]);
        ui.sliceSlider->setMaximum(imageDisplays.GetReferenceSize(m_currentSliceDir));
    }
}

void AIRWindow::ChangeInteractionMode() {
    if (ui.actionDrawing->isChecked()) {
        ui.graphicsView->drawingMode();
    } else {
        ui.graphicsView->transformMode();
    }
}