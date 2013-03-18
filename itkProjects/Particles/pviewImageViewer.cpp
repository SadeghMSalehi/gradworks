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


using namespace std;
using namespace pi;

class vtkMouseHandler: public vtkCommand {
private:
    ImageViewer* m_imageViewer;
public:
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
        viewMat->Print(cout);

        switch (eventId) {
            case vtkCommand::LeftButtonReleaseEvent: {
                m_imageViewer->m_sliceImg.SetVTKTransform(viewMat);
                m_imageViewer->ResampleSlice();
                style->OnLeftButtonUp();
                cout << "Left button up!" << endl;
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
        style->Print(cout);
        style->AddObserver(vtkCommand::LeftButtonReleaseEvent, m_mouseHandler);
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
        m_sliceImg.SetImage(fixedImg);
        RealImage::SizeType sz = fixedImg->GetBufferedRegion().GetSize();
        m_sliceImg.SelectSlice(2, sz[2] / 2);
        m_fixedPixmap = m_scene.addPixmap(m_sliceImg.GetPixmap());
        ui.fixedSliceSlider->setMaximum(sz[2]);
        ui.fixedSliceSlider->setValue(sz[2]/2);
    }
}

void ImageViewer::ResampleSlice() {
    if (m_fixedPixmap == NULL) {
        return;
    }
    m_sliceImg.GetResampled();
    m_fixedPixmap = m_scene.addPixmap(m_sliceImg.GetPixmap());
}

void ImageViewer::on_fixedOpacity_sliderMoved(int n) {
    if (m_fixedPixmap == NULL) {
        return;
    }
    m_fixedPixmap->setOpacity(ui.fixedOpacity->value()/255.0);
    m_fixedPixmap->update();
}

void ImageViewer::on_fixedSliceSlider_sliderMoved(int n) {
    m_sliceImg.SelectSlice(2, n);
    if (m_fixedPixmap != NULL) {
        m_scene.removeItem(m_fixedPixmap);
    }
    m_fixedPixmap = m_scene.addPixmap(m_sliceImg.GetPixmap());
}


void ImageViewer::on_zoomSlider_sliderMoved(int n) {
    QTransform transform;
    transform.scale(n/100.0, n/100.0);
    ui.graphicsView->setTransform(transform);
}
