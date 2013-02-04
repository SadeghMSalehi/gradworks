#include "dialogBSplineVis.h"
#include "iostream"
#include "QShowEvent"
#include "QCloseEvent"
#include "QMouseEvent"
#include "QGraphicsItem"

#include "myBSplineRegistration.h"
#include "vnlCommon.h"
#include "itkImageIO.h"
#include "itkCheckerBoardImageFilter.h"
#include "itkResampleImageFilter.h"
#include "QGraphicsGridItem.h"
#include "myImageTransform.h"
#include "myImageParticlesAlgorithm.h"
#include "mainwindow.h"

#include "QPrinter"
#include "QPainter"
#include "QFileDialog"

using namespace std;

BSplineVisDialog::BSplineVisDialog(QWidget* parent) : QDialog(parent) {
    ui.setupUi(this);
    QObject::connect(ui.showNoImage, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.showWarpedSlice, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.showCheckerboard, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.showWarpedCheckerboard, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.imageOpacity, SIGNAL(sliderMoved(int)), this, SLOT(updateScene()));
    QObject::connect(ui.showNoGrid, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.showCoordinateGrid, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.showDisplacementField, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.showDetJacobian, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.groupBox, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.showLandmarks, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.showWarpedLandmarks, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.landmarkSize, SIGNAL(valueChanged(double)), this, SLOT(updateScene()));
    QObject::connect(ui.gridResolution, SIGNAL(valueChanged(int)), this, SLOT(updateScene()));
    QObject::connect(ui.showLandmarks, SIGNAL(toggled(bool)), this, SLOT(updateScene()));

    ui.bspView->setScene(&m_Scene);

    m_Key = 0;
    m_Index = -1;
    m_Algo = NULL;

    double spacing[] = { 1, 1 };
    itkcmds::itkImageIO<SliceType> io;
    m_RefImage = io.NewImageT(100, 100, 1);
    m_RefImage->SetSpacing(spacing);
    
    ui.toolBox->setCurrentWidget(ui.landmarkControlPage);

    CreateGridAndCheckerboards(m_RefImage);
    updateScene();
}

BSplineVisDialog::~BSplineVisDialog() {

}


// construct grid and checkerboard image to compare before and after warping
// the coordinate system is defined by reference image
//
void BSplineVisDialog::CreateGridAndCheckerboards(SliceType::Pointer refImage) {
    itkcmds::itkImageIO<SliceType> io;

    m_RefImage = refImage;
    m_DetJacobian = io.NewImageT(refImage);

    VNLVector patterns(2);
    patterns[0] = patterns[1] = 10;
    m_SrcImage = ImageContainer::CreateCheckerBoards(m_RefImage, patterns);
    m_Field = DisplacementFieldType::Pointer(NULL);

    SliceType::SizeType refSize = refImage->GetBufferedRegion().GetSize();
    gX.set_size(refSize[0], refSize[1]);
    gY.set_size(refSize[0], refSize[1]);

    SliceType::IndexType idx;
    SliceType::PointType point;
    
    // assume gX and gY are physical coordinate system
    for (int i = 0; i < gX.rows(); i++) {
        for (int j = 0; j < gY.cols(); j++) {
            idx[0] = i;
            idx[1] = j;
            m_RefImage->TransformIndexToPhysicalPoint(idx, point);
            gX[i][j] = point[0];
            gY[i][j] = point[1];
        }
    }
}


void BSplineVisDialog::updateScene() {
    renderScene(m_Scene, false);
}

// refresh ui elements
//
void BSplineVisDialog::renderScene(QGraphicsScene& scene, bool isPrinter) {
    scene.clear();

    if (ui.showCheckerboard->isChecked()) {
        RGBAImageType::Pointer image = ImageContainer::CreateBitmap(m_SrcImage, ui.imageOpacity->value());
        QPixmap qPixmap = ImageContainer::CreatePixmap(image);
        scene.addPixmap(qPixmap)->setZValue(-10);
    }

    if (ui.showWarpedCheckerboard->isChecked()) {
        RGBAImageType::Pointer image = ImageContainer::CreateBitmap(m_DstImage, ui.imageOpacity->value());
        QPixmap qPixmap = ImageContainer::CreatePixmap(image);
        scene.addPixmap(qPixmap)->setZValue(-10);
    }

    if (ui.showWarpedSlice->isChecked()) {
        RGBAImageType::Pointer image = ImageContainer::CreateBitmap(m_WarpedSlice, ui.imageOpacity->value());
        QPixmap qPixmap = ImageContainer::CreatePixmap(image);
        scene.addPixmap(qPixmap)->setZValue(-10);

    }

    if (m_RefImage.IsNotNull()) {
        const int gridRes = ui.gridResolution->value();
        QPen blackPen = QPen(QColor::fromRgbF(0,0,0,0.5));
        blackPen.setWidthF(.1);
        QPen linePen(QColor::fromRgbF(.5,.5,.5,ui.imageOpacity->value() / 255.0));
        if (ui.showCoordinateGrid->isChecked()) {
            QGraphicsGridItem* originalGrid = new QGraphicsGridItem();
            if (isPrinter) {
                // can printer have alpha value?
                originalGrid->SetPen(blackPen);
            } else {
                originalGrid->SetPen(linePen);
            }
            originalGrid->SetResolution(gridRes);
            originalGrid->SetGrid(gX, gY);
            scene.addItem(originalGrid);
            ui.bspView->fitInView(originalGrid, Qt::KeepAspectRatio);
        }
        if (ui.showWarpedCoordinateGrid->isChecked()) {
            QGraphicsGridItem* warpedGrid = new QGraphicsGridItem();
            if (isPrinter) {
                // can printer have alpha value?
                if (!ui.showNoImage->isChecked()) {
                    blackPen.setColor(QColor::fromRgbF(1, 1, 0));
                    blackPen.setWidthF(0.2);
                }
                warpedGrid->SetPen(blackPen);
            } else {
                warpedGrid->SetPen(linePen);
            }
            warpedGrid->SetResolution(gridRes);
            warpedGrid->SetGrid(tX, tY);
            scene.addItem(warpedGrid);
            ui.bspView->fitInView(warpedGrid, Qt::KeepAspectRatio);
        }
    }


    if (ui.showDetJacobian->isChecked() && m_DetJacobian.IsNotNull()) {
        RGBAImageType::Pointer image = ImageContainer::CreateBitmap(m_DetJacobian);
        QPixmap qDetPixmap = ImageContainer::CreatePixmap(image);
        scene.addPixmap(qDetPixmap);
    }


    if (m_Field.IsNotNull() && ui.showDisplacementField->isChecked()) {
        FieldIteratorType iter(m_Field, m_Field->GetBufferedRegion());
        for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
            DisplacementFieldType::IndexType idx = iter.GetIndex();
            DisplacementFieldType::PointType point;
            m_Field->TransformIndexToPhysicalPoint(idx, point);
            VectorType v = iter.Get();
            if (v.GetNorm() < 1) {
                continue;
            }
            DisplacementFieldType::PointType point2;
            point2[0] = point[0] + v[0];
            point2[1] = point[1] + v[1];
            scene.addLine(point[0], point[1], point2[0], point2[1], QPen(Qt::yellow));
        }

    }



    if (ui.showWarpedLandmarks->isChecked()) {
        double sz = ui.landmarkSize->value();
        for (int i = 0; i < m_WarpedLandmarks.rows(); i++) {
            double x = m_WarpedLandmarks[i][0];
            double y = m_WarpedLandmarks[i][1];
            QGraphicsItem* p = scene.addRect(x - sz, y - sz, sz*2, sz*2, QPen(Qt::yellow), QBrush(Qt::yellow, Qt::SolidPattern));
        }
    }


    // drawing order is important for correct pick up
    // only show when showLandmarks is checked
    //
    if (ui.showLandmarks->isChecked()) {
        double sz = ui.landmarkSize->value();
        for (int i = 0; i < m_VectorList.size(); i++) {
            QRectF& xy = m_VectorList[i];
            QPen linePen(Qt::yellow);
            QPen sourcePen(Qt::red);
            QBrush sourceBrush(Qt::red, Qt::SolidPattern);
            QPen targetPen(Qt::blue);
            QBrush targetBrush(Qt::blue, Qt::SolidPattern);
            linePen.setWidthF(sz * .5);
            if (isPrinter) {
                linePen.setWidthF(sz);
                linePen.setColor(QColor::fromRgbF(0.3, 0.6, 0.3));
            }
            scene.addLine(xy.left(), xy.top(), xy.right(), xy.bottom(), linePen);
            QGraphicsItem* dx = scene.addEllipse(xy.left() - sz, xy.top() - sz, sz*2, sz*2, sourcePen, sourceBrush);
            QGraphicsItem* dy = scene.addEllipse(xy.right() - sz, xy.bottom() - sz, sz*2, sz*2, targetPen, targetBrush);

            dx->setData(1, i);
            dy->setData(2, i);
            
        }
    }


}

void BSplineVisDialog::on_bspView_mouseReleased(QMouseEvent *event) {
}

void BSplineVisDialog::on_bspView_mouseMoved(QMouseEvent *event) {
    if (event->buttons() == Qt::LeftButton) {
        QPointF pos = ui.bspView->mapToScene(event->pos());
        if (m_Key > 0) {
            switch (m_Key) {
                case 1:
                    m_VectorList[m_Index].setTopLeft(pos);
                    break;
                case 2:
                    m_VectorList[m_Index].setBottomRight(pos);
                    break;
                default:
                    break;
            }
            updateScene();
        }
    }
}


void BSplineVisDialog::on_bspView_mousePressed(QMouseEvent *event) {
    QPointF pos = ui.bspView->mapToScene(event->pos());
    QGraphicsItem* item = m_Scene.itemAt(pos);
    if (item != NULL) {
        if (!item->data(1).isNull()) {
            int key = item->data(1).toInt();
            m_Index = key;
            m_Key = 1;
        } else if (!item->data(2).isNull()) {
            int key = item->data(2).toInt();
            m_Index = key;
            m_Key = 2;
        } else {
            m_Index = -1;
            m_Key = 0;
        }
    } else {
        m_Index = -1;
        m_Key = 0;
    }
}


//
//void BSplineVisDialog::on_bspViewZoom_sliderMoved(int val) {
////    double zoom = ui.bspViewZoom->value() / 10.0f;
////    cout << "Zoom: " << zoom << endl;
//
//    QTransform transform = QTransform::fromScale(zoom, zoom);
//
//    ui.bspView->setTransform(transform);
//    updateScene();
//}

void BSplineVisDialog::on_addPairButton_clicked() {
    QRectF rect;
    rect.setCoords(10, 10, 20, 20);
    m_VectorList.push_back(rect);

    updateScene();
}

void BSplineVisDialog::on_copyPointsButton_clicked() {
    m_Algo = ((MainWindow*) parent())->GetImageParticlesAlgorithm();

    // select the current slice of the first image as reference image
    ImageContainer::Pointer img = m_Algo->GetImage(0);
    m_RefImage = img->GetSlice();

    // reset grid and checkerboard
    CreateGridAndCheckerboards(m_RefImage);

    // clear current landmarks
    m_VectorList.clear();

    const OptimizerParametersType& particles = m_Algo->GetCurrentParams();
    const int nPoints = m_Algo->GetNumberOfPoints();
    const int nParams = m_Algo->GetNumberOfParams();

    // copy landmarks from the first and second
    for (int i = 0; i < nPoints; i++) {
        QRectF rect;
        rect.setCoords(particles[SDim*i], particles[SDim*i+1], particles[nParams+SDim*i], particles[nParams+SDim*i+1]);
        m_VectorList.push_back(rect);
    }
    
    updateScene();
}

void BSplineVisDialog::on_updateField_clicked() {
//    m_Algo = ((MainWindow*) parent())->GetImageParticlesAlgorithm();

    VNLMatrix data(2, m_VectorList.size() * 2);
    for (int i = 0; i < m_VectorList.size(); i++) {
        data[0][i*2] = m_VectorList[i].left();
        data[0][i*2+1] = m_VectorList[i].top();
        data[1][i*2] = m_VectorList[i].right();
        data[1][i*2+1] = m_VectorList[i].bottom();
    }

    m_WarpedLandmarks.set_size(0, 0);

    if (ui.txfBspline->isChecked()) {
        my::BSplineRegistration breg;
        PropertyAccess props(this);
        breg.SetPropertyAccess(props);
        breg.SetLandmarks(m_VectorList.size(), data[0], data[1]);
        breg.SetReferenceImage(m_RefImage);
        breg.Update();
        FieldTransformType::Pointer txf = breg.GetTransform();

        m_Field = breg.GetDisplacementField();
        m_DetJacobian = breg.GetDeterminantOfJacobian();

        ImageContainer::WarpGrid(txf.GetPointer(), gX, gY, tX, tY);
        m_DstImage = ImageContainer::TransformSlice(m_SrcImage, txf.GetPointer());
        m_WarpedSlice = ImageContainer::TransformSlice(m_Algo->GetImage(1)->GetSlice(), txf.GetPointer());

        m_WarpedLandmarks.set_size(m_VectorList.size(), 2);

        FieldTransformType::InputPointType p;
        FieldTransformType::OutputPointType q;
        for (int i = 0; i < m_VectorList.size(); i++) {
            p[0] = m_VectorList[i].x();
            p[1] = m_VectorList[i].y();
            q = txf->TransformPoint(p);
            m_WarpedLandmarks[i][0] = q[0];
            m_WarpedLandmarks[i][1] = q[1];
        }

//        m_WarpedLabelSlice = ImageContainer::TransformSlice(m_Algo->GetImage(1)->GetLabelSliceAsSliceType(), txf.GetPointer(), true);        
    } else if (ui.txfBSplineFFD->isChecked()) {
        my::BSplineRegistration breg;
        PropertyAccess props(this);
        breg.SetPropertyAccess(props);
        breg.SetLandmarks(m_VectorList.size(), data[0], data[1]);
        breg.SetReferenceImage(m_RefImage);
        breg.SetUseFreeFormDeformation(true);
        breg.UpdateDeformation();
        SliceTransformType::Pointer txf = breg.GetFreeFormTransform();

        ImageContainer::WarpGrid(txf.GetPointer(), gX, gY, tX, tY);
//        cout << tX << ", " << tY << endl;
        m_DstImage = ImageContainer::TransformSlice(m_SrcImage, txf.GetPointer());
        m_WarpedSlice = ImageContainer::TransformSlice(m_Algo->GetImage(1)->GetSlice(), txf.GetPointer());
//        m_WarpedLabelSlice = ImageContainer::TransformSlice(m_Algo->GetImage(1)->GetLabelSliceAsSliceType(), txf.GetPointer(), true);
    } else {
        my::KernelTransformPointer transform;
        int type = 0;
        if (ui.txfTPS->isChecked()) {
            type = 0;
        } else if (ui.txfEBS->isChecked()) {
            type = 1;
        } else if (ui.txfTPSR2logR->isChecked()) {
            type = 2;
        }
        my::ImageTransform imageTxf;
        transform = imageTxf.CreateKernelTransform(type, m_VectorList.size(), data[0], data[1]);

        ImageContainer::WarpGrid(transform.GetPointer(), gX, gY, tX, tY);
        m_DstImage = ImageContainer::TransformSlice(m_SrcImage, transform.GetPointer());
        m_WarpedSlice = ImageContainer::TransformSlice(m_Algo->GetImage(1)->GetSlice(), transform.GetPointer());
//        m_WarpedLabelSlice = ImageContainer::TransformSlice(m_Algo->GetImage(1)->GetLabelSliceAsSliceType(), transform.GetPointer(), true);
    }

    ui.showWarpedCoordinateGrid->setChecked(true);
    updateScene();
}

void BSplineVisDialog::on_saveWarpedImage_clicked() {
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save NRRD File"), ".", tr("NRRD File (*.nrrd)"));

    if (fileName.isNull()) {
        return;
    }

    itkcmds::itkImageIO<SliceType> io;
    io.WriteImageT(fileName.toUtf8().data(), m_WarpedSlice);
}

void BSplineVisDialog::showEvent(QShowEvent* event) {

}

void BSplineVisDialog::closeEvent(QCloseEvent* event) {
    
}

void BSplineVisDialog::on_printButton_clicked() {
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save PDF File"), "/tmpfs", tr("PDF File (*.pdf)"));
    
    if (fileName.isNull()) {
        return;
    }

    QGraphicsScene printerScene;
    renderScene(printerScene, true);
    
    QPrinter printer(QPrinter::HighResolution);
    QSizeF paperSize;
    paperSize.setWidth(printerScene.sceneRect().width());
    paperSize.setHeight(printerScene.sceneRect().height());
    printer.setPaperSize(paperSize, QPrinter::Millimeter);
    printer.setPageMargins(5, 5, 5, 5, QPrinter::Millimeter);
    printer.setOrientation(QPrinter::Portrait);
    printer.setOutputFormat(QPrinter::PdfFormat);
    printer.setOutputFileName(fileName);


    QPainter painter(&printer);
    painter.setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing | QPainter::SmoothPixmapTransform | QPainter::NonCosmeticDefaultPen);
    printerScene.render(&painter);
    painter.end();
}
