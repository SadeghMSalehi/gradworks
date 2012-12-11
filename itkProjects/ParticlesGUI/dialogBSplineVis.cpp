#include "dialogBSplineVis.h"
#include "iostream"
#include "QShowEvent"
#include "QCloseEvent"
#include "QMouseEvent"
#include "QGraphicsItem"

#include "myBSplineRegistration.h"
#include "vnlCommon.h"
#include "itkImageIO.h"

using namespace std;

BSplineVisDialog::BSplineVisDialog(QWidget* parent) : QDialog(parent) {
    ui.setupUi(this);

    ui.bspView->setScene(&m_Scene);

    m_Key = 0;
    m_Index = -1;
    
    itkcmds::itkImageIO<SliceType> io;
    m_RefImage = io.NewImageT(100, 100, 1);
    
    updateScene();
}

BSplineVisDialog::~BSplineVisDialog() {

}

void BSplineVisDialog::updateScene() {
    m_Scene.clear();

    if (m_RefImage.IsNotNull()) {
        SliceType::SizeType sz = m_RefImage->GetBufferedRegion().GetSize();
        QPen pen(QColor::fromRgbF(1, 1, 1, 0.5));
        for (int i = 0; i <= sz[0]; i += 5) {
            QGraphicsItem* l = m_Scene.addLine(i, 0, i, sz[1], pen);
            l->setEnabled(false);
        }
        for (int j = 0; j <= sz[1]; j += 5) {
            QGraphicsItem* l = m_Scene.addLine(0, j, sz[0], j, pen);
            l->setEnabled(false);
        }
    }
    
    
    // drawing order is important for correct pick up
    double sz = 2 / 2;
    for (int i = 0; i < m_VectorList.size(); i++) {
        QRectF& xy = m_VectorList[i];
        m_Scene.addLine(xy.left(), xy.top(), xy.right(), xy.bottom(), QPen(Qt::white));
        QGraphicsItem* dx = m_Scene.addEllipse(xy.left() - sz, xy.top() - sz, sz*2, sz*2, QPen(Qt::red), QBrush(Qt::red, Qt::SolidPattern));
        dx->setData(1, i);
        QGraphicsItem* dy = m_Scene.addEllipse(xy.right() - sz, xy.bottom() - sz, sz*2, sz*2, QPen(Qt::blue), QBrush(Qt::blue, Qt::SolidPattern));
        dy->setData(2, i);

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

    cout << "Click position: " << pos.x() << ", " << pos.y() << endl;
}

void BSplineVisDialog::on_bspViewZoom_sliderMoved(int val) {
    double zoom = ui.bspViewZoom->value() / 10.0f;
    cout << "Zoom: " << zoom << endl;

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
    my::BSplineRegistration breg;

    VNLMatrix data(2, m_VectorList.size() * 2);
    for (int i = 0; i < m_VectorList.size(); i++) {
        data[0][i*2] = m_VectorList[i].left();
        data[0][i*2+1] = m_VectorList[i].top();
        data[1][i*2] = m_VectorList[i].right();
        data[1][i*2+1] = m_VectorList[i].bottom();
    }
    breg.SetLandmarks(m_VectorList.size(), data[0], data[1]);
    breg.SetReferenceImage(m_RefImage);
    
}

void BSplineVisDialog::showEvent(QShowEvent* event) {

}

void BSplineVisDialog::closeEvent(QCloseEvent* event) {
    
}


