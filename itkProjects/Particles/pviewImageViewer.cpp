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

QGraphicsVolumeItem::QGraphicsVolumeItem(QGraphicsItem* parent): QGraphicsPixmapItem(parent) {
    checkerBoardPattern = 0;
}

void QGraphicsVolumeItem::doPostProcess() {
    if (checkerBoardPattern == 0) {
        return;
    }
    
    QPixmap pmap = pixmap();
    QPainter p(&pmap);
    
    int w = pmap.width();
    int h = pmap.height();
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            
        }
    }
}

void QGraphicsVolumeItem::updatePixmap() {
    if (_resampledImg.IsNull()) {
        return;
    }

    typedef itk::ScalarToARGBColormapImageFilter<RealImage, RGBAVolumeType> ScalarToRGBFilter;
    typename ScalarToRGBFilter::Pointer rgbFilter = ScalarToRGBFilter::New();
    rgbFilter->SetInput(_resampledImg);
    rgbFilter->UseManualScalingOn();
    rgbFilter->UseIntensityWindowOn();
    rgbFilter->SetMinimumValue(_windowMin);
    rgbFilter->SetMaximumValue(_windowMax);
    rgbFilter->SetAlphaValue(255);
    rgbFilter->Update();

    RGBAVolumeType::Pointer rgbImg = rgbFilter->GetOutput();
    QImage qImg((unsigned char*) rgbImg->GetBufferPointer(), w, h, QImage::Format_ARGB32);
    setPixmap(QPixmap::fromImage(qImg));
    
    doPostProcess();
}

void QGraphicsVolumeItem::setImage(RealImage::Pointer img) {
    _srcImg = img;

    TransformType::ParametersType params, center;
    params.SetSize(12);
    center.SetSize(3);
    params.Fill(0);

    RealImage::RegionType gridRegion = _srcImg->GetBufferedRegion();
    IntIndex idx1 = gridRegion.GetIndex();
    IntIndex idx2 = gridRegion.GetUpperIndex();
    IntIndex centerIdx;
    for (int i = 0; i < 3; i++) {
        centerIdx[i] = (idx1[i]+idx2[i])/2.0;
    }
    RealImage::PointType centerPoint;
    _srcImg->TransformIndexToPhysicalPoint(centerIdx, centerPoint);
    for (int i = 0; i < 3; i++) {
        center[i] = centerPoint[i];
        // to make identity transform
        params[3*i+i] = 1;
    }

    _transform = TransformType::New();
    _transform->SetParameters(params);
    _transform->SetFixedParameters(center);
}

void QGraphicsVolumeItem::setResampleGrid(RealImage::Pointer grid) {
    _resampleGrid = grid;

    RealImage::RegionType region = _resampleGrid->GetBufferedRegion();
    // x-y, y-z, z-x

    w = region.GetSize(0);
    h = region.GetSize(1);

    if (_srcImg.IsNotNull()) {
        resampleGrid();
        updatePixmap();
    }
}

void QGraphicsVolumeItem::setTransform(vtkMatrix4x4 *mat) {
    if (_resampleGrid.IsNull() || _transform.IsNull()) {
        return;
    }

    TransformType::ParametersType params = _transform->GetParameters();
    params.Fill(0);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            params[3 * i + j] = mat->GetElement(i, j);
        }
    }
    _transform->SetParameters(params);

    resampleGrid();
    updatePixmap();
}

void QGraphicsVolumeItem::resampleGrid() {
    if (_resampleGrid.IsNull() || _srcImg.IsNull()) {
        return;
    }
    typedef itk::ResampleImageFilter<RealImage, RealImage> ResampleFilter;
    typename ResampleFilter::Pointer resampleFilter = ResampleFilter::New();
    typename ResampleFilter::TransformType* transformInput =
    dynamic_cast<typename ResampleFilter::TransformType*>(_transform.GetPointer());
    resampleFilter->SetInput(_srcImg);
    resampleFilter->SetReferenceImage(_resampleGrid);
    resampleFilter->UseReferenceImageOn();
    if (_transform.IsNotNull()) {
        resampleFilter->SetTransform(transformInput);
    }
    resampleFilter->Update();
    _resampledImg = resampleFilter->GetOutput();

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
                    m_imageViewer->m_movingItem->setTransform(viewMat);
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
    m_fixedItem = m_movingItem = NULL;

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

    m_fixedItem = new QGraphicsVolumeItem();
    m_movingItem = new QGraphicsVolumeItem();

    renderer->ResetCamera();
    renderer->Render();
}

ImageViewer::~ImageViewer() {
    delete m_mouseHandler;
}

void ImageViewer::LoadImage(QString fileName) {
    if (!fileName.isEmpty()) {
        ImageIO<RealImage> io;
        fixedImg = io.ReadImage(fileName.toStdString());

        typedef itk::StatisticsImageFilter<RealImage> StatFilter;
        StatFilter::Pointer statFilter = StatFilter::New();
        statFilter->SetInput(fixedImg);
        statFilter->Update();

        DataReal pixMin = statFilter->GetMinimum();
        DataReal pixMax = statFilter->GetMaximum();

        ui.intensitySlider->setRealMin(pixMin);
        ui.intensitySlider->setRealMax(pixMax);

        ui.intensitySlider->setLowValue(ui.intensitySlider->minimum());
        ui.intensitySlider->setHighValue(ui.intensitySlider->maximum());
        
        m_fixedItem->setImage(fixedImg);
        m_fixedItem->setWindowRange(pixMin, pixMax);
        m_scene.addItem(m_fixedItem);

        SetFixedSlice(2);

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

void ImageViewer::LoadMovingImage(QString fileName) {
    if (!fileName.isEmpty()) {
        ImageIO<RealImage> io;
        movingImg = io.ReadImage(fileName.toStdString());

        if (movingImg.IsNull()) {
            return;
        }

        typedef itk::StatisticsImageFilter<RealImage> StatFilter;
        StatFilter::Pointer statFilter = StatFilter::New();
        statFilter->SetInput(movingImg);
        statFilter->Update();

        DataReal pixMin = statFilter->GetMinimum();
        DataReal pixMax = statFilter->GetMaximum();

        ui.intensitySlider2->setEnabled(true);
        ui.intensitySlider2->setRealMin(pixMin);
        ui.intensitySlider2->setRealMax(pixMax);

        ui.intensitySlider2->setLowValue(ui.intensitySlider2->minimum());
        ui.intensitySlider2->setHighValue(ui.intensitySlider2->maximum());

        m_movingItem->setImage(movingImg);
        m_movingItem->setWindowRange(pixMin, pixMax);
        m_scene.addItem(m_movingItem);

        SetMovingSlice(2);
    }
}

void ImageViewer::SetFixedSlice(int dir, int idx) {
    if (fixedImg.IsNull()) {
        return;
    }
    RealImage::SizeType sz = fixedImg->GetBufferedRegion().GetSize();
    ui.fixedSliceSlider->setMaximum(sz[2]);
    if (idx == -1) {
        idx = sz[2] / 2.0;
    }
    ui.fixedSliceSlider->setValue(idx);

    RealImage::RegionType resampleRegion = fixedImg->GetBufferedRegion();
    RealImage::IndexType idx1 = resampleRegion.GetIndex();
    RealImage::IndexType idx2 = resampleRegion.GetUpperIndex();
    idx1[dir] = idx;
    idx2[dir] = idx;
    resampleRegion.SetIndex(idx1);
    resampleRegion.SetUpperIndex(idx2);

    typedef itk::ExtractImageFilter<RealImage, RealImage> ExtractFilterType;
    typename ExtractFilterType::Pointer filter = ExtractFilterType::New();
    filter->SetInput(fixedImg);
    filter->SetExtractionRegion(resampleRegion);
    filter->Update();
    fixedGrid = filter->GetOutput();

    m_fixedItem->setResampleGrid(fixedGrid);
}

void ImageViewer::SetMovingSlice(int dir, int idx) {
    if (movingImg.IsNull()) {
        return;
    }
    m_movingItem->setResampleGrid(fixedGrid);
}

void ImageViewer::on_fixedOpacity_sliderMoved(int n) {
    if (m_fixedItem == NULL) {
        return;
    }
    m_fixedItem->setOpacity(ui.fixedOpacity->value()/255.0);
    m_fixedItem->update();
}

void ImageViewer::on_movingOpacity_sliderMoved(int n) {
    if (m_movingItem == NULL) {
        return;
    }
    m_movingItem->setOpacity(ui.movingOpacity->value()/255.0);
    m_movingItem->update();
}


void ImageViewer::on_fixedSliceSlider_sliderMoved(int n) {
    SetFixedSlice(2, n);
    SetMovingSlice(2);
}

void ImageViewer::on_zoomSlider_sliderMoved(int n) {
    QTransform transform;
    transform.scale(n/100.0, n/100.0);
    ui.graphicsView->setTransform(transform);
}

void ImageViewer::on_intensitySlider_lowValueChanged(int n) {
    m_fixedItem->setWindowRange(ui.intensitySlider->realLowValue(), ui.intensitySlider->realHighValue());
}

void ImageViewer::on_intensitySlider_highValueChanged(int n) {
    m_fixedItem->setWindowRange(ui.intensitySlider->realLowValue(), ui.intensitySlider->realHighValue());
}

void ImageViewer::on_intensitySlider2_lowValueChanged(int n) {
    m_movingItem->setWindowRange(ui.intensitySlider2->realLowValue(), ui.intensitySlider2->realHighValue());
}

void ImageViewer::on_intensitySlider2_highValueChanged(int n) {
    m_movingItem->setWindowRange(ui.intensitySlider2->realLowValue(), ui.intensitySlider2->realHighValue());
}

void ImageViewer::on_fixedMask_checked(bool check) {
    RealImage::Pointer img = m_fixedItem->getImage();

}

void ImageViewer::on_movingMask_checked(bool check) {
    
}

