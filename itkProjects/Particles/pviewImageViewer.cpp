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


using namespace std;
using namespace pi;

typedef itk::Image<RGBAPixel, 3> RGBAVolumeType;

void QGraphicsImageItem::updatePixmap() {
    if (_srcImg.IsNull()) {
        return;
    }

    typedef itk::ScalarToARGBColormapImageFilter<RealImage, RGBAVolumeType> ScalarToRGBFilter;
    typename ScalarToRGBFilter::Pointer rgbFilter = ScalarToRGBFilter::New();
    rgbFilter->SetInput(_srcImg);
    rgbFilter->UseManualScalingOn();
    rgbFilter->UseIntensityWindowOn();
    rgbFilter->SetMinimumValue(_windowMin);
    rgbFilter->SetMaximumValue(_windowMax);
    rgbFilter->SetAlphaValue(255);
    rgbFilter->Update();

    RGBAVolumeType::Pointer rgbImg = rgbFilter->GetOutput();
    QImage qImg((unsigned char*) rgbImg->GetBufferPointer(), w, h, QImage::Format_ARGB32);
    setPixmap(QPixmap::fromImage(qImg));
}

void QGraphicsImageItem::setImage(RealImage::Pointer img) {
    _srcImg = img;
    // XXX: have to change if the slice is not an axial.
    RealImage::SizeType sz = _srcImg->GetBufferedRegion().GetSize();
    w = sz[0];
    h = sz[1];
    updatePixmap();
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
                    m_imageViewer->m_movingImg.SetVTKTransform(viewMat);
                    m_imageViewer->ResampleSlice();
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
    m_fixedPixmap = m_movingPixmap = NULL;

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
}

ImageViewer::~ImageViewer() {
    delete m_mouseHandler;
}

void ImageViewer::LoadImage(QString fileName) {
    if (!fileName.isEmpty()) {
        ImageIO<RealImage> io;
        fixedImg = io.ReadImage(fileName.toStdString());
        m_fixedImg.SetImage(fixedImg);
        RealImage::SizeType sz = fixedImg->GetBufferedRegion().GetSize();
        m_fixedImg.SetSliceAsGrid(2, sz[2]/2);
        m_fixedPixmap = m_scene.addPixmap(m_fixedImg.GetPixmap());
        ui.fixedSliceSlider->setMaximum(sz[2]);
        ui.fixedSliceSlider->setValue(sz[2]/2);

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
    }
}

void ImageViewer::LoadMovingImage(QString fileName) {
    if (!fileName.isEmpty()) {
        ImageIO<RealImage> io;
        RealImage::Pointer img = io.ReadImage(fileName.toStdString());
        m_movingImg.SetImage(img);
        m_movingImg.SetResampleGrid(m_fixedImg.GetResampleGrid());
        m_movingPixmap = m_scene.addPixmap(m_movingImg.GetPixmap());
    }
}

void ImageViewer::ResampleSlice() {
    if (m_movingImg.HasNoImage()) {
        return;
    }
    if (m_movingPixmap != NULL) {
        m_scene.removeItem(m_movingPixmap);
    }
    m_movingImg.Resample();
    m_movingPixmap = m_scene.addPixmap(m_movingImg.GetPixmap());
    m_movingPixmap->setOpacity(ui.movingOpacity->value()/255.0);
    m_movingPixmap->update();
}

void ImageViewer::on_fixedOpacity_sliderMoved(int n) {
    if (m_fixedPixmap == NULL) {
        return;
    }
    m_fixedPixmap->setOpacity(ui.fixedOpacity->value()/255.0);
    m_fixedPixmap->update();
}

void ImageViewer::on_movingOpacity_sliderMoved(int n) {
    if (m_movingPixmap == NULL) {
        return;
    }
    m_movingPixmap->setOpacity(ui.movingOpacity->value()/255.0);
    m_movingPixmap->update();
}


void ImageViewer::on_fixedSliceSlider_sliderMoved(int n) {
    m_fixedImg.SetSliceAsGrid(2, n);
    if (m_fixedPixmap != NULL) {
        m_scene.removeItem(m_fixedPixmap);
    }
    m_fixedPixmap = m_scene.addPixmap(m_fixedImg.GetPixmap());
    m_fixedPixmap->setOpacity(ui.fixedOpacity->value()/255.0);
    m_fixedPixmap->update();
}


void ImageViewer::on_zoomSlider_sliderMoved(int n) {
    QTransform transform;
    transform.scale(n/100.0, n/100.0);
    ui.graphicsView->setTransform(transform);
}

void ImageViewer::on_intensitySlider_lowValueChanged(int n) {
    cout << "Low: " << ui.intensitySlider->realLowValue() << endl;
}

void ImageViewer::on_intensitySlider_highValueChanged(int n) {
    cout << "High: " << ui.intensitySlider->realHighValue() << endl;
}
