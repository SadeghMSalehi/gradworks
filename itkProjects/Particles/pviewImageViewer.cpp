//
//  pviewImageTransform.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/14/13.
//
//

#include "pviewImageViewer.h"
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
#include "QGraphicsItem"


using namespace std;
using namespace pi;

typedef itk::Image<RGBAPixel, 3> RGBAVolumeType;

bool QGraphicsCompositeImageItem::CheckCompositeBuffer(int id1) {
    if (_imageDisplays == NULL) {
        return false;
    }
    
    pi::ImageIO<RealImage> io;
    ImageDisplayType* img1 = _imageDisplays->at(id1);
    
    // check foreground image
    RealImage::Pointer fImg = img1->GetResampled(_resampleIdx);
    if (fImg.IsNull()) {
        return false;
    }
    
    // check composite image
    if (_compositeImage.IsNull()
        || _compositeImage->GetBufferedRegion() != fImg->GetBufferedRegion()) {
        _compositeImage = io.CopyImage(fImg);
    }
    return _compositeImage.IsNotNull();
}

void QGraphicsCompositeImageItem::CompositeAlpha(int id1, int id2) {
    if (!CheckCompositeBuffer(id1)) {
        return;
    }
    
    // prepare input and output buffer
    ImageDisplayType* img1 = _imageDisplays->at(id1);
    RealImage::Pointer fImg = img1->GetResampled(_resampleIdx);
    const int nElems = fImg->GetPixelContainer()->Size();
    DataReal* cBuf = _compositeImage->GetBufferPointer();
    DataReal* fBuf = fImg->GetBufferPointer();
    
    if (id2 >= 0) {
        // if there is more than one image,
        ImageDisplayType* img2 = _imageDisplays->at(id2);
        RealImage::Pointer mImg = _imageDisplays->at(_movingId)->GetResampled(_resampleIdx);
        DataReal* mBuf = mImg->GetBufferPointer();
        for (int i = 0; i < nElems; i++) {
            const DataReal f = Clamp(fBuf[i], img1->histogram.rangeMin, img1->histogram.rangeMax);
            const DataReal m = Clamp(mBuf[i], img2->histogram.rangeMin, img2->histogram.rangeMax);
            cBuf[i] = _alpha*f+(1-_alpha)*m;
        }
    } else {
        // for the case of single image
        for (int i = 0; i < nElems; i++) {
            const DataReal f = Clamp(fBuf[i], img1->histogram.rangeMin, img1->histogram.rangeMax);
            cBuf[i] = _alpha * f;
        }
    }
}

void QGraphicsCompositeImageItem::CompositeCheckerBoard(int id1, int id2) {
    if (!CheckCompositeBuffer(id1)) {
        return;
    }
    
    if (id2 < 0 || id2 >= _imageDisplays->Count()) {
        return;
    }
    
    // prepare input and output buffer
    ImageDisplayType* img1 = _imageDisplays->at(id1);
    RealImage::Pointer fImg = img1->GetResampled(_resampleIdx);
    DataReal* fBuf = fImg->GetBufferPointer();

    ImageDisplayType* img2 = _imageDisplays->at(id2);
    RealImage::Pointer mImg = _imageDisplays->at(_movingId)->GetResampled(_resampleIdx);
    DataReal* mBuf = mImg->GetBufferPointer();
    
    DataReal* cBuf = _compositeImage->GetBufferPointer();
    
    int w = fImg->GetBufferedRegion().GetSize(0);
    int h = fImg->GetBufferedRegion().GetSize(1);
    
    const int cbszW = w / _cbCols;
    const int cbszH = h / _cbRows;
    int k = 0;
    for (int i = 0; i < w; i++) {
        for (int j = 0; j < h; j++, k++) {
            const DataReal f = Clamp(fBuf[k], img1->histogram.rangeMin, img1->histogram.rangeMax);
            const DataReal m = Clamp(mBuf[k], img2->histogram.rangeMin, img2->histogram.rangeMax);
            if ((i/cbszW)%2 == (j/cbszH)%2) {
                cBuf[k] = f;
            } else {
                cBuf[k] = m;
            }
        }
    }
}

void QGraphicsCompositeImageItem::Refresh(int id1, int id2) {
    if (_imageDisplays == NULL || _imageDisplays->Count() <= id1 || id1 < 0) {
        return;
    }

    if (_compositionMode == Alpha) {
        CompositeAlpha(id1, id2);
    } else if (_compositionMode == CheckerBoard) {
        CompositeCheckerBoard(id1, id2);
    }

    // if the composition image is not generated,
    if (_compositeImage.IsNull()) {
        return;
    }

    // convert to color image
    // assume that the composite image ranges from 0 to 65535 (in short range)
    ColorFilterType::Pointer colorFilter = ColorFilterType::New();
    colorFilter->SetInput(_compositeImage);
    colorFilter->UseManualScalingOn();
    colorFilter->SetMinimumValue(0);
    colorFilter->SetMaximumValue(65535);
    colorFilter->Update();
    _rgbImage = colorFilter->GetOutput();
    
    // may have to correct to deal with various slice
    int w = _rgbImage->GetBufferedRegion().GetSize(0);
    int h = _rgbImage->GetBufferedRegion().GetSize(1);

    setPixmap(QPixmap::fromImage(QImage((unsigned char*) _rgbImage->GetBufferPointer(), w, h, QImage::Format_ARGB32)));
    update();
}

class vtkMouseHandler: public vtkCommand {
private:
    ImageViewer* m_imageViewer;
    bool m_moving;
public:
    vtkMouseHandler() {
        m_moving = false;
    }

    void SetImageViewer(ImageViewer* imageViewer) {
        m_imageViewer = imageViewer;
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
                    if (m_imageViewer->imageDisplays.Count() > 1) {
                        m_imageViewer->UpdateMovingDisplayTransform(viewMat);
                        m_imageViewer->UpdateCompositeDisplay();
                    }
                }
                style->OnMouseMove();
                break;
            }
            default:
                return;
        }
    }
};


ImageViewer::ImageViewer(QWidget* parent) {
    ui.setupUi(this);
    vtkRenderer* renderer = vtkRenderer::New();
    m_mouseHandler = new vtkMouseHandler();
    m_mouseHandler->SetImageViewer(this);

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
}

ImageViewer::~ImageViewer() {
    delete m_mouseHandler;
}

void ImageViewer::LoadImage(QString fileName) {
    if (!fileName.isEmpty()) {
        ImageIO<RealImage> io;
        RealImage::Pointer img = io.ReadImage(fileName.toStdString());

        imageDisplays.AddImage(img);
        if (imageDisplays.Count() == 1) {
            m_fixedId = 0;
            OnFixedImageLoaded();
        } else {
            m_movingId = 1;
            OnMovingImageLoaded();
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

void ImageViewer::OnFixedImageLoaded() {
    ImageDisplayType& dispImg = imageDisplays[m_fixedId];
    imageDisplays.SetReferenceId(m_fixedId);
    
    // slice axis selection
    ui.fixedSliceSlider->setMaximum(dispImg.srcImg->GetBufferedRegion().GetSize()[2]);
    ui.fixedSliceSlider->setValue(dispImg.srcImg->GetBufferedRegion().GetSize()[2]/2.0);
    
    ui.intensitySlider->setRealMin(dispImg.histogram.dataMin);
    ui.intensitySlider->setRealMax(dispImg.histogram.dataMax);
    
    ui.intensitySlider->setLowValue(ui.intensitySlider->minimum());
    ui.intensitySlider->setHighValue(ui.intensitySlider->maximum());
    
    imageDisplays.SetSliceGrid(2, ui.fixedSliceSlider->value());
}

void ImageViewer::OnMovingImageLoaded() {
    ImageDisplayType& dispImg = imageDisplays[m_movingId];
    
    ui.intensitySlider2->setRealMin(dispImg.histogram.dataMin);
    ui.intensitySlider2->setRealMax(dispImg.histogram.dataMax);
    
    ui.intensitySlider2->setLowValue(ui.intensitySlider2->minimum());
    ui.intensitySlider2->setHighValue(ui.intensitySlider2->maximum());
    
    ui.intensitySlider2->setEnabled(true);
}

void ImageViewer::UpdateMovingDisplayTransform(vtkMatrix4x4* mat) {
    if (m_movingId >= 0 && imageDisplays.Count() > m_movingId) {
        imageDisplays[m_movingId].SetAffineTransform(mat);
    }
}

void ImageViewer::UpdateCompositeDisplay() {
    m_compositeDisplay->Refresh(m_fixedId, m_movingId);
}

void ImageViewer::on_compositeOpacity_sliderMoved(int n) {
    double alpha = 1-double(n)/ui.compositeOpacity->maximum();
    m_compositeDisplay->CompositionModeToAlpha(alpha);
    UpdateCompositeDisplay();
}


void ImageViewer::on_fixedSliceSlider_sliderMoved(int n) {
    imageDisplays.SetSliceGrid(2, n);
    UpdateCompositeDisplay();
}

void ImageViewer::on_zoomSlider_sliderMoved(int n) {
    QTransform transform;
    transform.scale(n/100.0, n/100.0);
    ui.graphicsView->setTransform(transform);
}

void ImageViewer::on_intensitySlider_lowValueChanged(int n) {
    on_intensitySlider_highValueChanged(0);
}

void ImageViewer::on_intensitySlider_highValueChanged(int n) {
    // correct to handle multiple image histogram intensity
    imageDisplays[m_fixedId].histogram.rangeMin = ui.intensitySlider->realLowValue();
    imageDisplays[m_fixedId].histogram.rangeMax = ui.intensitySlider->realHighValue();
    UpdateCompositeDisplay();
}

void ImageViewer::on_intensitySlider2_lowValueChanged(int n) {
    on_intensitySlider2_highValueChanged(0);
}

void ImageViewer::on_intensitySlider2_highValueChanged(int n) {
    // correct to handle multiple image histogram intensity
    imageDisplays[m_movingId].histogram.rangeMin = ui.intensitySlider2->realLowValue();
    imageDisplays[m_movingId].histogram.rangeMax = ui.intensitySlider2->realHighValue();
    UpdateCompositeDisplay();
}

void ImageViewer::on_alphaOptions_toggled(bool b) {
    ui.checkerBoardOptions->setChecked(false);
    on_compositeOpacity_sliderMoved(ui.compositeOpacity->value());
}

void ImageViewer::on_checkerBoardOptions_toggled(bool b) {
    ui.alphaOptions->setChecked(false);
    m_compositeDisplay->CompositionModeToCheckerBoard(ui.checkerBoardRows->value(), ui.checkerBoardCols->value());
    UpdateCompositeDisplay();
}