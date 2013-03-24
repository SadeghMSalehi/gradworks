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

void QGraphicsCompositeImageItem::CompositeAlpha() {
    if (_imageDisplays == NULL) {
        return;
    }

    pi::ImageIO<RealImage> io;
    RealImage::Pointer fImg = _imageDisplays->at(_fixedId)->GetResampled(_resampleIdx);
    if (fImg.IsNull()) {
        cout << "no resampled image" << endl;
        return;
    }
    if (_compositeImage.IsNull()
        || _compositeImage->GetBufferedRegion() != fImg->GetBufferedRegion()) {
        _compositeImage = io.CopyImage(fImg);
    }
    int nElems = fImg->GetPixelContainer()->Size();

    DataReal* fBuf = fImg->GetBufferPointer();
    DataReal* cBuf = _compositeImage->GetBufferPointer();

    
    if (_imageDisplays->Count() <= _movingId || _imageDisplays->at(_movingId)->srcImg.IsNull()) {
        for (int i = 0; i < nElems; i++) {
            const DataReal f = Clamp(fBuf[i], _imageDisplays->at(_fixedId)->histogram.rangeMin, _imageDisplays->at(_fixedId)->histogram.rangeMax);
            cBuf[i] = _alpha * f;
        }
    } else {
        RealImage::Pointer mImg = _imageDisplays->at(_movingId)->GetResampled(_resampleIdx);
        DataReal* mBuf = mImg->GetBufferPointer();
        for (int i = 0; i < nElems; i++) {
            const DataReal f = Clamp(fBuf[i], _imageDisplays->at(_fixedId)->histogram.rangeMin, _imageDisplays->at(_fixedId)->histogram.rangeMax);
            const DataReal m = Clamp(mBuf[i], _imageDisplays->at(_movingId)->histogram.rangeMin, _imageDisplays->at(_movingId)->histogram.rangeMax);
            cBuf[i] = _alpha*f+(1-_alpha)*m;
        }
    }
}

void QGraphicsCompositeImageItem::Refresh() {
    if (_imageDisplays == NULL || _imageDisplays->Count() <= _fixedId || _fixedId < 0) {
        return;
    }

    if ((_movingId < 0 || _movingId >= _imageDisplays->Count()) && (_compositionMode == CheckerBoard)) {
        return;
    }

    if (_compositionMode == Alpha) {
        CompositeAlpha();
    } else if (_compositionMode == CheckerBoard) {

    }

    if (_compositeImage.IsNull()) {
        return;
    }

    ColorFilterType::Pointer colorFilter = ColorFilterType::New();
    colorFilter->SetInput(_compositeImage);
    colorFilter->UseManualScalingOn();
    colorFilter->SetMinimumValue(0);
    colorFilter->SetMaximumValue(65535);
    colorFilter->Update();

    _rgbImage = colorFilter->GetOutput();
    int w = _rgbImage->GetBufferedRegion().GetSize(0);
    int h = _rgbImage->GetBufferedRegion().GetSize(1);

    setPixmap(QPixmap::fromImage(QImage((unsigned char*) _rgbImage->GetBufferPointer(), w, h, QImage::Format_ARGB32)));
    update();
}

//QGraphicsVolumeItem::QGraphicsVolumeItem(QGraphicsItem* parent): QGraphicsPixmapItem(parent) {
//}
//
//void QGraphicsVolumeItem::doPostProcess() {
//    if (true) {
//        return;
//    }
//    
//    QPixmap pmap = pixmap();
//    QPainter p(&pmap);
//    
//    int w = pmap.width();
//    int h = pmap.height();
//
//    QBrush brush(QColor::fromRgb(0, 0, 0), Qt::SolidPattern);
//    int r = 4, c = 4;
//    for (int i = 0; i < r; i++) {
//        for (int j = 0; j < c; j++) {
//            if (i%2 == j%2) {
//                int rw = w/c;
//                if (j*(w/c)+w/c >= w) {
//                    rw = w - j*(w/c);
//                }
//                int rh = h/r;
//                if (i*(h/r)+h/r >= h) {
//                    rh = h - i*(h/r);
//                }
//                p.eraseRect(j*(w/c), i*(h/r), rw, rh);
//                p.fillRect(j*(w/c), i*(h/r), rw, rh, brush);
//            }
//        }
//    }
//
//    setPixmap(pmap);
//}

//void QGraphicsVolumeItem::updatePixmap() {
//    if (_resampledImg.IsNull()) {
//        return;
//    }
//
//    typedef itk::ScalarToARGBColormapImageFilter<RealImage, RGBAVolumeType> ScalarToRGBFilter;
//    typename ScalarToRGBFilter::Pointer rgbFilter = ScalarToRGBFilter::New();
//    rgbFilter->SetInput(_resampledImg);
//    rgbFilter->UseManualScalingOn();
//    rgbFilter->UseIntensityWindowOn();
//    rgbFilter->SetMinimumValue(_windowMin);
//    rgbFilter->SetMaximumValue(_windowMax);
//    rgbFilter->SetAlphaValue(255);
//    rgbFilter->Update();
//
//    RGBAVolumeType::Pointer rgbImg = rgbFilter->GetOutput();
//    QImage qImg((unsigned char*) rgbImg->GetBufferPointer(), w, h, QImage::Format_ARGB32);
//    setPixmap(QPixmap::fromImage(qImg));
//    
//    doPostProcess();
//}
//
//void QGraphicsVolumeItem::setImage(RealImage::Pointer img) {
//    _srcImg = img;
//
//    TransformType::ParametersType params, center;
//    params.SetSize(12);
//    center.SetSize(3);
//    params.Fill(0);
//
//    RealImage::RegionType gridRegion = _srcImg->GetBufferedRegion();
//    IntIndex idx1 = gridRegion.GetIndex();
//    IntIndex idx2 = gridRegion.GetUpperIndex();
//    IntIndex centerIdx;
//    for (int i = 0; i < 3; i++) {
//        centerIdx[i] = (idx1[i]+idx2[i])/2.0;
//    }
//    RealImage::PointType centerPoint;
//    _srcImg->TransformIndexToPhysicalPoint(centerIdx, centerPoint);
//    for (int i = 0; i < 3; i++) {
//        center[i] = centerPoint[i];
//        // to make identity transform
//        params[3*i+i] = 1;
//    }
//
//    _transform = TransformType::New();
//    _transform->SetParameters(params);
//    _transform->SetFixedParameters(center);
//}
//
//void QGraphicsVolumeItem::setResampleGrid(RealImage::Pointer grid) {
//    _resampleGrid = grid;
//
//    RealImage::RegionType region = _resampleGrid->GetBufferedRegion();
//    // x-y, y-z, z-x
//
//    w = region.GetSize(0);
//    h = region.GetSize(1);
//
//    if (_srcImg.IsNotNull()) {
//        resampleGrid();
//        updatePixmap();
//    }
//}
//
//void QGraphicsVolumeItem::setTransform(vtkMatrix4x4 *mat) {
//    if (_resampleGrid.IsNull() || _transform.IsNull()) {
//        return;
//    }
//
//    TransformType::ParametersType params = _transform->GetParameters();
//    params.Fill(0);
//    for (int i = 0; i < 3; i++) {
//        for (int j = 0; j < 3; j++) {
//            params[3 * i + j] = mat->GetElement(i, j);
//        }
//    }
//    _transform->SetParameters(params);
//
//    resampleGrid();
//    updatePixmap();
//}
//
//void QGraphicsVolumeItem::resampleGrid() {
//    if (_resampleGrid.IsNull() || _srcImg.IsNull()) {
//        return;
//    }
//    typedef itk::ResampleImageFilter<RealImage, RealImage> ResampleFilter;
//    typename ResampleFilter::Pointer resampleFilter = ResampleFilter::New();
//    typename ResampleFilter::TransformType* transformInput =
//    dynamic_cast<typename ResampleFilter::TransformType*>(_transform.GetPointer());
//    resampleFilter->SetInput(_srcImg);
//    resampleFilter->SetReferenceImage(_resampleGrid);
//    resampleFilter->UseReferenceImageOn();
//    if (_transform.IsNotNull()) {
//        resampleFilter->SetTransform(transformInput);
//    }
//    resampleFilter->Update();
//    _resampledImg = resampleFilter->GetOutput();
//
//}
//
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
//                    m_imageViewer->m_movingItem->setTransform(viewMat);
                    viewMat->Print(cout);
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

        ImageDisplayType& dispImg = imageDisplays.AddImage(img);
        
        if (imageDisplays.Count() == 1) {
            imageDisplays.SetReferenceId(0);

            // slice axis selection
            ui.fixedSliceSlider->setMaximum(dispImg.srcImg->GetBufferedRegion().GetSize()[2]);
            ui.fixedSliceSlider->setValue(dispImg.srcImg->GetBufferedRegion().GetSize()[2]/2.0);

            ui.intensitySlider->setRealMin(dispImg.histogram.dataMin);
            ui.intensitySlider->setRealMax(dispImg.histogram.dataMax);

            ui.intensitySlider->setLowValue(ui.intensitySlider->minimum());
            ui.intensitySlider->setHighValue(ui.intensitySlider->maximum());

            imageDisplays.SetSliceGrid(2, ui.fixedSliceSlider->value());
        }
        m_compositeDisplay->Refresh();

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

void ImageViewer::on_compositeOpacity_sliderMoved(int n) {
    double alpha = 1-double(n)/ui.compositeOpacity->maximum();
    m_compositeDisplay->CompositionModeToAlpha(alpha);
    m_compositeDisplay->Refresh();
}


void ImageViewer::on_fixedSliceSlider_sliderMoved(int n) {
    imageDisplays.SetSliceGrid(2, n);
    m_compositeDisplay->Refresh();    
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
    imageDisplays[0].histogram.rangeMin = ui.intensitySlider->realLowValue();
    imageDisplays[0].histogram.rangeMax = ui.intensitySlider->realHighValue();
    m_compositeDisplay->Refresh();
}

