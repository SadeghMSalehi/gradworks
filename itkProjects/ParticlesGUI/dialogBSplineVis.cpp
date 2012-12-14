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

using namespace std;

BSplineVisDialog::BSplineVisDialog(QWidget* parent) : QDialog(parent) {
    ui.setupUi(this);

    QObject::connect(ui.showOriginal, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.imageOpacity, SIGNAL(sliderMoved(int)), this, SLOT(updateScene()));
    QObject::connect(ui.showCoordinateGrid, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.showDisplacementField, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.showDetJacobian, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.groupBox, SIGNAL(toggled(bool)), this, SLOT(updateScene()));

    ui.bspView->setScene(&m_Scene);

    m_Key = 0;
    m_Index = -1;

    double spacing[] = { 1, 1 };
    itkcmds::itkImageIO<SliceType> io;
    m_RefImage = io.NewImageT(100, 100, 1);
    m_RefImage->SetSpacing(spacing);

    m_DetJacobian = io.NewImageT(100, 100, 1);
    m_DetJacobian->SetSpacing(spacing);
    VNLVector patterns(2);
    patterns[0] = patterns[1] = 50;
    m_SrcImage = ImageContainer::CreateCheckerBoards(m_RefImage, patterns);
    m_Field = DisplacementFieldType::Pointer(NULL);

    gX.set_size(100, 100);
    gY.set_size(100, 100);

    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < 100; j++) {
            gX[i][j] = i;
            gY[i][j] = j;
        }
    }

    updateScene();
}

BSplineVisDialog::~BSplineVisDialog() {

}

void BSplineVisDialog::updateScene() {
    m_Scene.clear();

    if (ui.showOriginal->isChecked()) {
        RGBAImageType::Pointer image = ImageContainer::CreateBitmap(m_SrcImage, ui.imageOpacity->value());
        QPixmap qPixmap = ImageContainer::CreatePixmap(image);
        m_Scene.addPixmap(qPixmap)->setZValue(-10);
    }

    if (ui.showTransformed->isChecked()) {
        RGBAImageType::Pointer image = ImageContainer::CreateBitmap(m_DstImage, ui.imageOpacity->value());
        QPixmap qPixmap = ImageContainer::CreatePixmap(image);
        m_Scene.addPixmap(qPixmap)->setZValue(-10);
    }

    if (m_RefImage.IsNotNull()) {
        const int gridRes = ui.gridResolution->value();
        if (ui.showCoordinateGrid->isChecked()) {
            QPen pen(QColor::fromRgbF(1, 1, 1, ui.imageOpacity->value() / 255.0));
            QGraphicsGridItem* originalGrid = new QGraphicsGridItem();
            originalGrid->SetPen(pen);
            originalGrid->SetResolution(gridRes);
            originalGrid->SetGrid(gX, gY);
            m_Scene.addItem(originalGrid);
        }
        if (ui.showWarpedCoordinateGrid->isChecked() && m_Field.IsNotNull()) {
            QPen pen(QColor::fromRgbF(1, 1, 1, ui.imageOpacity->value() / 255.0));
            QGraphicsGridItem* warpedGrid = new QGraphicsGridItem();
            warpedGrid->SetPen(pen);
            warpedGrid->SetResolution(gridRes);
            warpedGrid->SetGrid(tX, tY);
            m_Scene.addItem(warpedGrid);
        }
    }


    if (ui.showDetJacobian->isChecked() && m_DetJacobian.IsNotNull()) {
        RGBAImageType::Pointer image = ImageContainer::CreateBitmap(m_DetJacobian);
        QPixmap qDetPixmap = ImageContainer::CreatePixmap(image);
        m_Scene.addPixmap(qDetPixmap);
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
            m_Scene.addLine(point[0], point[1], point2[0], point2[1], QPen(Qt::yellow));
        }

    }


    // drawing order is important for correct pick up
    double sz = 1.5;
    for (int i = 0; i < m_VectorList.size(); i++) {
        QRectF& xy = m_VectorList[i];
        m_Scene.addLine(xy.left(), xy.top(), xy.right(), xy.bottom(), QPen(Qt::white));
        QGraphicsItem* dx = m_Scene.addEllipse(xy.left() - sz, xy.top() - sz, sz*2, sz*2, QPen(Qt::red), QBrush(Qt::red, Qt::SolidPattern));
        dx->setData(1, i);
        QGraphicsItem* dy = m_Scene.addEllipse(xy.right() - sz, xy.bottom() - sz, sz*2, sz*2, QPen(Qt::blue), QBrush(Qt::blue, Qt::SolidPattern));
        dy->setData(2, i);

    }

    QRectF bounds(0, 0, 100, 100);
    ui.bspView->fitInView(bounds, Qt::KeepAspectRatio);
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

    cout << "Click position: " << pos.x() << ", " << pos.y() << endl;
}

void BSplineVisDialog::on_bspViewZoom_sliderMoved(int val) {
    double zoom = ui.bspViewZoom->value() / 10.0f;
//    cout << "Zoom: " << zoom << endl;

    QTransform transform = QTransform::fromScale(zoom, zoom);

    ui.bspView->setTransform(transform);
    updateScene();
}

void BSplineVisDialog::on_addPairButton_clicked() {
    QRectF rect;
    rect.setCoords(10, 10, 20, 20);
    m_VectorList.push_back(rect);

    updateScene();
}

void BSplineVisDialog::on_updateField_clicked() {
    VNLMatrix data(2, m_VectorList.size() * 2);
    for (int i = 0; i < m_VectorList.size(); i++) {
        data[0][i*2] = m_VectorList[i].left();
        data[0][i*2+1] = m_VectorList[i].top();
        data[1][i*2] = m_VectorList[i].right();
        data[1][i*2+1] = m_VectorList[i].bottom();
    }
    
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
    }
    updateScene();
}

void BSplineVisDialog::showEvent(QShowEvent* event) {

}

void BSplineVisDialog::closeEvent(QCloseEvent* event) {
    
}


