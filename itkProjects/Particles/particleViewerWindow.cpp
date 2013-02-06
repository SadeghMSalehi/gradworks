#include "particleViewerWindow.h"
#include "myParticleCore.h"
#include "myParticleBSpline.h"
#include "QFileDialog"
#include "QGraphicsGridItem.h"

pi::ParticleSystem g_System;
pi::VNLMatrix g_GX, g_GY, g_TX, g_TY;

ParticleViewerWindow::ParticleViewerWindow(QWidget* parent) {
    ui.setupUi(this);
    ui.graphicsView->setScene(&m_Scene);

    QObject::connect(ui.actionShow_Coordinate_Grid, SIGNAL(triggered(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.actionShow_Particles, SIGNAL(triggered(bool)), this, SLOT(updateScene()));

}

ParticleViewerWindow::~ParticleViewerWindow() {

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
    g_System.LoadSystem(file);
    createGrid();
    updateScene();
}

void ParticleViewerWindow::updateScene() {
    m_Scene.clear();

    const pi::ParticleSubjectArray& subjects = g_System.GetSubjects();
    const int nSubjects = subjects.size();
    const int nPoints = subjects[0].GetNumberOfPoints();

    if (ui.actionShow_Coordinate_Grid->isChecked()) {
        QGraphicsGridItem* warpedGrid = new QGraphicsGridItem();
        QPen linePen(QColor::fromRgbF(.5,.5,.5,128/255.0));
        warpedGrid->SetPen(linePen);
        warpedGrid->SetResolution(1);
        warpedGrid->SetGrid(g_TX, g_TY);
        m_Scene.addItem(warpedGrid);
    }

    if (ui.actionShow_Particles->isChecked()) {
        for (int n = 0; n < nSubjects; n++) {
            if (n > 0) {
                for (int i = 0; i < nPoints; i++) {
                    const pi::Particle& src = subjects[0][i];
                    const pi::Particle& dst = subjects[n][i];
                    m_Scene.addLine(src.x[0], src.x[1], dst.x[0], dst.x[1], QPen(Qt::green));
                }
            }
            for (int i = 0; i < nPoints; i++) {
                const pi::Particle& par = subjects[n][i];
                double x = par.x[0];
                double y = par.x[1];

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
