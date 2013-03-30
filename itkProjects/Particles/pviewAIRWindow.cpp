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
#include "qutils.h"
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QMimeData>
#include <QList>
#include <QUrl>
#include <QShortcut>
#include "airAlgorithmManager.h"

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
    ui.progressBar->hide();

    // set up for segmentation tools
    _drawingToolButtons = new QButtonGroup(this);
    _drawingToolButtons->addButton(ui.freePath);
    _drawingToolButtons->addButton(ui.freeBrush);
    _drawingToolButtons->addButton(ui.eraseBrush);
    ui.freePath->setChecked(true);


    ui.backgroundBrush->addItem("Any", QVariant(255));
    ui.backgroundBrush->addItem("Empty", QVariant(0));
    for (int i = 1; i < __maxLabelColors; i++) {
        QImage iconImg(16,16, QImage::Format_RGB32);
        iconImg.fill(ui.graphicsView->colorTable(uchar(i)));
        QIcon icon(QPixmap::fromImage(iconImg));
        ui.foregroundBrush->addItem(icon, QString("%1").arg(i), QVariant(i));
        ui.backgroundBrush->addItem(icon, QString("%1").arg(i), QVariant(i));
    }


    QList<QAction*> segmentationMenuActions;
    segmentationMenuActions.append(ui.actionLoadSegmentation);
    segmentationMenuActions.append(ui.actionSaveSegmentation);
    segmentationMenuActions.append(ui.actionNewSegmentation);
    segmentationMenuActions.append(ui.actionPropagateLabel);
    ui.segmentationMenu->addActions(segmentationMenuActions);

    ui.copyButton->setDefaultAction(ui.actionCopyLabel);
    ui.pasteButton->setDefaultAction(ui.actionPasteLabel);


    // transformation widget
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

    pivtk::PolyDataPointer sphere = pivtk::CreateSphere(60, 60);
    m_propScene.SetRenderer(renderer);
    m_propScene.AddPolyData("sphere", sphere);
    m_propScene.SetRepresentation(VTK_WIREFRAME);

    renderer->ResetCamera();
    renderer->Render();


    // slice direction and navigation
    ui.graphicsView->setScene(&m_scene);

    m_fixedId = -1;
    m_movingId = -1;
    m_currentSliceDir = IJ;
    m_requestedSliceDir = IJ;

    setupShortcutKeys();


    m_compositeDisplay = new QGraphicsCompositeImageItem();
    m_compositeDisplay->SetImageDisplays(&imageDisplays);

    // Slice Direction
    m_sliceDirectionActions.addAction(ui.actionSliceIJ);
    m_sliceDirectionActions.addAction(ui.actionSliceJK);
    m_sliceDirectionActions.addAction(ui.actionSliceKI);
    m_sliceDirectionActions.setExclusive(true);
    ui.actionSliceIJ->setChecked(true);


#pragma mark QObject Signal-Slot Connections
    // image IO
    connect(this, SIGNAL(fileDropped(QString&)), this, SLOT(on_image1Name_fileDropped(QString&)));

    // slice view selection related signals
    connect(m_compositeDisplay, SIGNAL(translationChanged()), this, SLOT(UpdateTranslationWidget()));
    connect(ui.actionSliceIJ, SIGNAL(triggered()), this, SLOT(ChangeSliceDirection()));
    connect(ui.actionSliceJK, SIGNAL(triggered()), this, SLOT(ChangeSliceDirection()));
    connect(ui.actionSliceKI, SIGNAL(triggered()), this, SLOT(ChangeSliceDirection()));
    connect(ui.actionDrawing, SIGNAL(triggered()), this, SLOT(ChangeInteractionMode()));

    // working-set related signals
    connect(ui.actionNewWorkingSet, SIGNAL(triggered()), ui.multipleSliceView, SLOT(createWorkingSet()));
    connect(ui.actionClearWorkingSet, SIGNAL(triggered()), ui.multipleSliceView, SLOT(clearWorkingSet()));

    // drawing related signals
    connect(ui.freeBrush, SIGNAL(toggled(bool)), this, SLOT(ChangeInteractionMode()));
    connect(ui.freePath, SIGNAL(toggled(bool)), this, SLOT(ChangeInteractionMode()));
    connect(ui.eraseBrush, SIGNAL(toggled(bool)), this, SLOT(ChangeInteractionMode()));

    connect(ui.actionLoadSegmentation, SIGNAL(triggered()), this, SLOT(LoadSegmentation()));
    connect(ui.actionSaveSegmentation, SIGNAL(triggered()), this, SLOT(SaveSegmentation()));
    connect(ui.actionNewSegmentation, SIGNAL(triggered()), ui.graphicsView, SLOT(segmentationCleared()));
    connect(ui.actionPropagateLabel, SIGNAL(triggered()), this, SLOT(PropagateSegmentation()));
    connect(ui.labelOpacity, SIGNAL(sliderMoved(int)), ui.graphicsView, SLOT(labelOpacityChanged(int)));
    connect(ui.labelOpacity, SIGNAL(valueChanged(int)), ui.graphicsView, SLOT(labelOpacityChanged(int)));
    connect(ui.foregroundBrush, SIGNAL(currentIndexChanged(int)), this, SLOT(brushLabelChanged(int)));
    connect(ui.backgroundBrush, SIGNAL(currentIndexChanged(int)), this, SLOT(brushLabelChanged(int)));
    connect(ui.actionCopyLabel, SIGNAL(triggered()), this, SLOT(copyLabel()));
    connect(ui.actionPasteLabel, SIGNAL(triggered()), ui.graphicsView, SLOT(pasteLabel()));

    // slice navigation
    connect(ui.sliceNumber, SIGNAL(valueChanged(int)), ui.sliceSlider, SLOT(setValue(int)));
    connect(ui.sliceSlider, SIGNAL(valueChanged(int)), ui.sliceNumber, SLOT(setValue(int)));
    connect(ui.sliceSlider, SIGNAL(valueChanged(int)), ui.multipleSliceView, SLOT(currentSliceChanged(int)));
    connect(ui.multipleSliceView, SIGNAL(sliceDoubleClicked(int)), ui.sliceSlider, SLOT(setValue(int)));

    ui.multipleSliceView->addAction(ui.actionNewWorkingSet);
    ui.multipleSliceView->addAction(ui.actionPropagateLabel);
    ui.multipleSliceView->addAction(ui.actionClearWorkingSet);

    ui.multipleSliceView->hide();

    // Prepare algorithm manager
    _algorithms = new air::AlgorithmManager(this, ui);

}

AIRWindow::~AIRWindow() {
    m_mouseHandler->UnRegister();
}

bool AIRWindow::IsImage1Loaded() {
    return imageDisplays.IsValidId(m_fixedId);
}

bool AIRWindow::IsImage2Loaded() {
    return imageDisplays.IsValidId(m_movingId);
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
        QString fileName = __fileManager.openFile(QFileManager::Single, this, "Load a fixed image");
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
        QString fileName = __fileManager.openFile(QFileManager::Single, this, "Load a moving image");
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

    ui.graphicsView->createLabelVolumeIfNecessary(imageDisplays[0].srcImg);

    __fileManager.putFile(QFileManager::Single, fileName);

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
    ChangeSliceDirection();

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

void AIRWindow::on_actionResample_triggered() {
    if (imageDisplays.IsValidId(m_movingId)) {
        QString fileName = __fileManager.saveFile(QFileManager::Single, this, "Save Resampled");
        if (fileName.isNull()) {
            return;
        }
        pi::ImageIO<AIRImage> io;
        AIRImage::Pointer resampledImg = imageDisplays[m_movingId].Resample3D(imageDisplays.GetReferenceGrid(), 0);
        io.WriteImage(fileName.toUtf8().data(), resampledImg);
    }
}

void AIRWindow::on_actionLoadTransform_triggered() {
    if (imageDisplays.IsValidId(m_movingId)) {
        QString fileName = __fileManager.openFile(QFileManager::Single, this, "Load a transform");
        if (fileName.isNull()) {
            return;
        }
        pi::ImageIO<AIRImage> io;
        pi::ImageIO<AIRImage>::TransformType::Pointer transform = io.ReadTransform(fileName.toUtf8().data());
        imageDisplays[m_movingId].SetAffineTransform(transform);
        UpdateTranslationWidget();
        UpdateCompositeDisplay();
    }
}

void AIRWindow::on_actionSaveTransform_triggered() {
    if (imageDisplays.IsValidId(m_movingId)) {
        QString fileName = __fileManager.saveFile(QFileManager::Single, this, "Save a transform");
        if (fileName.isNull()) {
            return;
        }

        itk::TransformFileWriter::Pointer writer = itk::TransformFileWriter::New();
        writer->SetInput(imageDisplays[m_movingId].GetAffineTransform());
        writer->SetFileName(fileName.toUtf8().data());
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
            ui.multipleSliceView->currentSliceChanged(ui.sliceSlider->value());
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
    if (!imageDisplays.IsValidId(m_fixedId)) {
        return;
    }
    ChangeSliceIndex(n);
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

void AIRWindow::ChangeSliceIndex(int n) {
    m_currentSliceIndex[m_currentSliceDir] = n;
    imageDisplays.SetReferenceSlice(m_currentSliceDir, m_currentSliceIndex[m_currentSliceDir]);
    ui.graphicsView->sliceChanged(m_currentSliceDir, m_currentSliceIndex[m_currentSliceDir]);
}

void AIRWindow::ChangeSliceDirection() {
    if (imageDisplays.IsValidId(m_fixedId)) {
        m_currentSliceDir = Unknown;
        if (ui.actionSliceIJ->isChecked()) {
            m_currentSliceDir = IJ;
        } else if (ui.actionSliceJK->isChecked()) {
            m_currentSliceDir = JK;
        } else if (ui.actionSliceKI->isChecked()) {
            m_currentSliceDir = KI;
        }

        const int sliderValue = ui.sliceSlider->value();
        const int sliderMaximum = ui.sliceSlider->maximum();
        const int requestedSliceIdx = m_currentSliceIndex[m_currentSliceDir];
        const int requestedSliceMaximum = imageDisplays.GetReferenceSize(m_currentSliceDir)-1;

        if (requestedSliceMaximum < sliderMaximum) {
            ui.sliceSlider->setValue(m_currentSliceIndex[m_currentSliceDir]);
        }
        if (requestedSliceIdx > sliderMaximum) {
            ui.sliceNumber->setMaximum(imageDisplays.GetReferenceSize(m_currentSliceDir)-1);
            ui.sliceSlider->setMaximum(imageDisplays.GetReferenceSize(m_currentSliceDir)-1);
        }
        
        ui.sliceSlider->setValue(m_currentSliceIndex[m_currentSliceDir]);
        ui.sliceNumber->setMaximum(imageDisplays.GetReferenceSize(m_currentSliceDir)-1);
        ui.sliceSlider->setMaximum(imageDisplays.GetReferenceSize(m_currentSliceDir)-1);
        if (requestedSliceIdx == sliderValue) {
            on_sliceSlider_valueChanged(requestedSliceIdx);
        }
    }
}

void AIRWindow::ChangeInteractionMode() {
    if (ui.actionDrawing->isChecked()) {
        ui.labelTools->setEnabled(true);
        if (ui.freePath->isChecked()) {
            ui.graphicsView->drawingMode(QGraphicsGuideView::FreePath);
        } else if (ui.freeBrush->isChecked()) {
            ui.graphicsView->drawingMode(QGraphicsGuideView::FreeBrush);
        } else if (ui.eraseBrush->isChecked()) {
            ui.graphicsView->drawingMode(QGraphicsGuideView::EraseBrush);
        }
        ui.toolBox->setCurrentWidget(ui.segmentationTools);
    } else {
        ui.labelTools->setEnabled(false);
        ui.graphicsView->transformMode();
    }
}


void AIRWindow::LoadSegmentation() {
    QString fileName = __fileManager.openFile(QFileManager::Single, this, "Load Segmentation");
    if (fileName.isEmpty()) {
        return;
    }
    if (imageDisplays.IsValidId(m_fixedId)) {
        ui.graphicsView->loadLabelVolume(fileName, imageDisplays[m_fixedId].srcImg);
    }
}


void AIRWindow::SaveSegmentation() {
    QString fileName = __fileManager.saveFile(QFileManager::Single, this, "Save Segmentation");
    if (fileName.isEmpty()) {
        return;
    }
    ui.graphicsView->saveLabelVolume(fileName);
}

void AIRWindow::PropagateSegmentation() {
    std::vector<int> workingSet = ui.multipleSliceView->getWorkingSet();
    ui.graphicsView->propagateLabel(workingSet);
}

void AIRWindow::brushLabelChanged(int n) {
    int f = ui.foregroundBrush->itemData(ui.foregroundBrush->currentIndex()).value<int>();
    int b = ui.backgroundBrush->itemData(ui.backgroundBrush->currentIndex()).value<int>();
    ui.graphicsView->setInteractionLabels(f,b);
}


#pragma mark -
#pragma mark Private Methods

void AIRWindow::previousSlice() {
    ui.sliceSlider->setValue(ui.sliceSlider->value()-1);
}

void AIRWindow::nextSlice() {
    ui.sliceSlider->setValue(ui.sliceSlider->value()+1);
}

void AIRWindow::sliceZoomIn() {
    ui.zoomSlider->setValue(ui.zoomSlider->value() + ui.zoomSlider->pageStep());
}

void AIRWindow::sliceZoomOut() {
    ui.zoomSlider->setValue(ui.zoomSlider->value() - ui.zoomSlider->pageStep());
}

void AIRWindow::labelOpacityUp() {
    ui.labelOpacity->setValue(ui.labelOpacity->value() + ui.labelOpacity->pageStep());
}

void AIRWindow::labelOpacityDown() {
    ui.labelOpacity->setValue(ui.labelOpacity->value() - ui.labelOpacity->pageStep());
}

void AIRWindow::copyLabel() {
    ui.statusbar->showMessage(QString("Label #%1 Copied to Clipboard").arg(ui.sliceSlider->value()), 5000);
    ui.graphicsView->copyLabel();
}

void AIRWindow::setupShortcutKeys() {
    // navigation
    QShortcut* prevSliceKey = new QShortcut(QKeySequence(tr("a")), this);
    QShortcut* nextSliceKey = new QShortcut(QKeySequence(tr("d")), this);
    prevSliceKey->setContext(Qt::ApplicationShortcut);
    nextSliceKey->setContext(Qt::ApplicationShortcut);

    // zoom in out
    QShortcut* zoomInKey = new QShortcut(QKeySequence(tr("w")), this);
    QShortcut* zoomOutKey = new QShortcut(QKeySequence(tr("s")), this);
    zoomInKey->setContext(Qt::ApplicationShortcut);
    zoomOutKey->setContext(Qt::ApplicationShortcut);


    // zoom in out
    QShortcut* labelOpacityUp = new QShortcut(QKeySequence(tr("e")), this);
    QShortcut* labelOpacityDown = new QShortcut(QKeySequence(tr("q")), this);
    labelOpacityUp->setContext(Qt::ApplicationShortcut);
    labelOpacityDown->setContext(Qt::ApplicationShortcut);

    QShortcut* selectPathTool = new QShortcut(QKeySequence(tr("1")), this);
    QShortcut* selectBrushTool = new QShortcut(QKeySequence(tr("2")), this);
    QShortcut* selectEraseTool = new QShortcut(QKeySequence(tr("3")), this);
    QShortcut* selectCopyTool = new QShortcut(QKeySequence(tr("c")), this);
    QShortcut* selectPasteTool = new QShortcut(QKeySequence(tr("v")), this);

    selectPathTool->setContext(Qt::ApplicationShortcut);
    selectBrushTool->setContext(Qt::ApplicationShortcut);
    selectEraseTool->setContext(Qt::ApplicationShortcut);
    selectCopyTool->setContext(Qt::ApplicationShortcut);
    selectPasteTool->setContext(Qt::ApplicationShortcut);

    connect(prevSliceKey, SIGNAL(activated()), this, SLOT(previousSlice()));
    connect(nextSliceKey, SIGNAL(activated()), this, SLOT(nextSlice()));
    connect(zoomInKey, SIGNAL(activated()), this, SLOT(sliceZoomIn()));
    connect(zoomOutKey, SIGNAL(activated()), this, SLOT(sliceZoomOut()));
    connect(labelOpacityUp, SIGNAL(activated()), this, SLOT(labelOpacityUp()));
    connect(labelOpacityDown, SIGNAL(activated()), this, SLOT(labelOpacityDown()));

    connect(selectPathTool, SIGNAL(activated()), ui.freePath, SLOT(toggle()));
    connect(selectBrushTool, SIGNAL(activated()), ui.freeBrush, SLOT(toggle()));
    connect(selectEraseTool, SIGNAL(activated()), ui.eraseBrush, SLOT(toggle()));
    connect(selectCopyTool, SIGNAL(activated()), this, SLOT(copyLabel()));
    connect(selectPasteTool, SIGNAL(activated()), ui.graphicsView, SLOT(pasteLabel()));
}