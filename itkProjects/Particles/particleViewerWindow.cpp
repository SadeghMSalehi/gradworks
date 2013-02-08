#include "particleViewerWindow.h"
#include "myParticleCore.h"
#include "myParticleBSpline.h"
#include "QFileDialog"
#include "QGraphicsGridItem.h"

namespace pi {
    
    class ParticleSlice {
    public:
        typedef std::vector<Particle*> ParticlePointerVector;

        ParticleSlice() {}
        ~ParticleSlice() {};

        const int GetNumberOfPointsInSlice(int slice, int subj);
        const ParticlePointerVector& Get(int slice, int subj);
        void Update(int sliceDim, ParticleSubjectArray& shapes, LabelImage::Pointer refImage);

    private:
        boost::numeric::ublas::matrix<ParticlePointerVector> m_ParticlePointerMatrix;
    };
    static ParticleSlice::ParticlePointerVector emptyPoint;



    void ParticleSlice::Update(int sliceDim, ParticleSubjectArray& subjects, LabelImage::Pointer labelImage) {
        LabelImage::SizeType sz = labelImage->GetBufferedRegion().GetSize();

        const int nSubj = subjects.size();
        const int nPoints = subjects[0].GetNumberOfPoints();

        m_ParticlePointerMatrix.resize(sz[sliceDim], nSubj);
        for (int i = 0; i < nSubj; i++) {
            for (int j = 0; j < nPoints; j++) {
                for (int k = 0; k < sz[sliceDim]; k++) {
                    if (subjects[i][j].x[sliceDim] >= k && subjects[i][j].x[sliceDim] < k + 1) {
                        m_ParticlePointerMatrix(k, i).push_back(&subjects[i][j]);
                        break;
                    }
                }
            }
        }
    }

    const int ParticleSlice::GetNumberOfPointsInSlice(int slice, int subj) {
        if (m_ParticlePointerMatrix.size1() > slice && m_ParticlePointerMatrix.size2() > subj) {
            return m_ParticlePointerMatrix(slice, subj).size();
        }
        return 0;
    }

    const ParticleSlice::ParticlePointerVector&  ParticleSlice::Get(int slice, int subj) {
        if (m_ParticlePointerMatrix.size1() > slice && m_ParticlePointerMatrix.size2() > subj) {
            return m_ParticlePointerMatrix(slice, subj);
        }
        return emptyPoint;
    }
}

pi::ParticleSystem g_System;
pi::ParticleSlice g_Slice;
pi::VNLMatrix g_GX, g_GY, g_TX, g_TY;

ParticleViewerWindow::ParticleViewerWindow(QWidget* parent) {
    ui.setupUi(this);
    ui.graphicsView->setScene(&m_Scene);

    m_CurrentDirection = 1;
    fordim (k) {
        m_CurrentSlice[k] = -1;
    }

    QObject::connect(ui.actionShow_Coordinate_Grid, SIGNAL(triggered(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.actionShow_Particles, SIGNAL(triggered(bool)), this, SLOT(updateScene()));

}

ParticleViewerWindow::~ParticleViewerWindow() {

}

void ParticleViewerWindow::setupSlice() {
    pi::LabelImage::Pointer label = g_System.GetImageContext().GetLabel(0);
    pi::LabelImage::SizeType sz = label->GetBufferedRegion().GetSize();
    fordim (k) {
        m_Size[k] = sz[k];
        m_CurrentSlice[k] = m_Size[k]/2.0;
    }
    ui.sliceSlider->setValue(m_CurrentSlice[m_CurrentDirection]);
    try {
        g_Slice.Update(m_CurrentDirection, g_System.GetSubjects(), label);
    } catch (std::exception& e) {
        cout << e.what() << endl;
    }
}

void ParticleViewerWindow::createGrid() {
    pi::LabelImage::Pointer refImage = g_System.GetImageContext().GetLabel(0);
    pi::LabelImage::SizeType refSize = refImage->GetBufferedRegion().GetSize();
    g_GX.set_size(refSize[0], refSize[1]);
    g_GY.set_size(refSize[0], refSize[1]);

    // assume gX and gY are physical coordinate system
    for (int i = 0; i < g_GX.rows(); i++) {
        for (int j = 0; j < g_GY.cols(); j++) {
            g_GX[i][j] = i;
            g_GY[i][j] = j;
        }
    }

    // create transform
    const pi::ParticleSubjectArray& subjects = g_System.GetSubjects();
    pi::ParticleBSpline bspline;
    bspline.SetReferenceImage(refImage);
    bspline.EstimateTransform(subjects[1], subjects[0]);
    pi::FieldTransformType::Pointer transform = bspline.GetTransform();

    // warp grid
    const int nRows = g_GX.rows();
    const int nCols = g_GX.cols();
    g_TX.set_size(nRows, nCols);
    g_TY.set_size(nRows, nCols);
    for (int i = 0; i < nRows; i++) {
        for (int j = 0; j < nCols; j++) {
            pi::FieldTransformType::InputPointType inPoint;
            inPoint[0] = g_GX[i][j];
            inPoint[1] = g_GY[i][j];
            pi::FieldTransformType::OutputPointType outPoint;
            outPoint = transform->TransformPoint(inPoint);
            g_TX[i][j] = outPoint[0];
            g_TY[i][j] = outPoint[1];
        }
    }
}

void ParticleViewerWindow::load(char* file) {
    if (g_System.LoadSystem(file)) {
        setupSlice();
        createGrid();
        try {
            updateScene();
        } catch (std::exception& e) {
            cout << e.what() << endl;
        }
    }
}

void ParticleViewerWindow::updateScene() {
    m_Scene.clear();

    const pi::ParticleSubjectArray& subjects = g_System.GetSubjects();
    const int currentSlice = m_CurrentSlice[m_CurrentDirection];
    const int dim1 = (m_CurrentDirection + 1) % __Dim;
    const int dim2 = (m_CurrentDirection + __Dim - 1) % __Dim;
    const int nSubjects = subjects.size();
    const int nPoints = g_Slice.GetNumberOfPointsInSlice(currentSlice, 0);

    if (ui.actionShow_Coordinate_Grid->isChecked()) {
        QGraphicsGridItem* warpedGrid = new QGraphicsGridItem();
        QPen linePen(QColor::fromRgbF(.5,.5,.5,128/255.0));
        warpedGrid->SetPen(linePen);
        warpedGrid->SetResolution(3);
        warpedGrid->SetGrid(g_TX, g_TY);
        m_Scene.addItem(warpedGrid);
    }

    if (ui.actionShow_Particles->isChecked()) {
        for (int n = 0; n < nSubjects; n++) {
            if (n > 0) {
                // draw correspondence line
                for (int i = 0; i < nPoints; i++) {
                    const pi::Particle* src = g_Slice.Get(currentSlice, 0)[i];
                    const int idx = src->idx;
                    const pi::Particle& dst = subjects[n][idx];
                    m_Scene.addLine(src->x[dim1], src->x[dim2], dst.x[dim1], dst.x[dim2], QPen(Qt::green));
                }
            }

            // draw particels
            for (int i = 0; i < nPoints; i++) {
                const pi::Particle* refPar = g_Slice.Get(currentSlice, 0)[i];
                const int idx = refPar->idx;
                const pi::Particle& par = subjects[n][idx];
                double x = par.x[dim1];
                double y = par.x[dim2];
                // better if we have random color palette
                QColor pointColor = Qt::black;
                if (n == 0) {
                    pointColor = Qt::red;
                } else if (n == 1) {
                    pointColor = Qt::blue;
                } else if (n == 2) {
                    pointColor = Qt::yellow;
                }
                m_Scene.addEllipse(x-.5, y-.5, 1, 1, QPen(pointColor), QBrush(pointColor, Qt::SolidPattern));
            }
        }
    }

    QRectF sceneRect = m_Scene.sceneRect();
    ui.graphicsView->fitInView(sceneRect, Qt::KeepAspectRatio);
}

void ParticleViewerWindow::on_action_Open_triggered() {
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), ".", tr("Text file (*.txt)"));

    if (fileName.isNull()) {
        return;
    }
    
    load(fileName.toUtf8().data());
}

void ParticleViewerWindow::on_action_Close_triggered() {
    exit(1);
}

void ParticleViewerWindow::on_sliceSlider_valueChanged(int slice) {
    m_CurrentSlice[m_CurrentDirection] = slice;
    try {
        updateScene();
    } catch (std::exception& e) {
        cout << e.what() << endl;
    }
}