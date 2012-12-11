#include "dialogBSplineVis.h"
#include "iostream"
#include "QShowEvent"
#include "QCloseEvent"
#include "QMouseEvent"
#include "QGraphicsItem"

using namespace std;

BSplineVisDialog::BSplineVisDialog(QWidget* parent) : QDialog(parent) {
    ui.setupUi(this);

    ui.bspView->setScene(&m_Scene);

    m_XY = NULL;
    m_Key = 0;
}

BSplineVisDialog::~BSplineVisDialog() {

}

void BSplineVisDialog::updateScene() {
    m_Scene.clear();

    QBrush brush(Qt::yellow, Qt::SolidPattern);
    for (int i = 0; i < m_VectorList.size(); i++) {
        QRectF& xy = m_VectorList[i];
        QGraphicsItem* dx = m_Scene.addEllipse(xy.left() - .5, xy.top() - .5, 1, 1, QPen(Qt::yellow), brush);
        dx->setData(1, i);
        QGraphicsItem* dy = m_Scene.addEllipse(xy.right() - .5, xy.bottom() - .5, 1, 1, QPen(Qt::yellow), brush);
        dy->setData(2, i);
        QGraphicsItem* dl = m_Scene.addLine(xy.left(), xy.top(), xy.right(), xy.bottom(), QPen(Qt::white));
        dl->setData(3, i);

    }
}

void BSplineVisDialog::on_bspView_mouseReleased(QMouseEvent *event) {
}

void BSplineVisDialog::on_bspView_mouseMoved(QMouseEvent *event) {
    if (event->buttons() == Qt::LeftButton) {
        QPointF pos = ui.bspView->mapToScene(event->pos());
    }
}


void BSplineVisDialog::on_bspView_mousePressed(QMouseEvent *event) {
    QPointF pos = ui.bspView->mapToScene(event->pos());
    QGraphicsItem* item = m_Scene.itemAt(pos);
    if (item != NULL) {
        if (!item->data(1).isNull()) {
            int key = item->data(1).toInt();
            m_XY = &m_VectorList[key];
            m_Key = 1;
        } else if (!item->data(2).isNull()) {
            int key = item->data(2).toInt();
            m_XY = &m_VectorList[key];
            m_Key = 2;
        } else if (!item->data(3).isNull()) {
            int key = item->data(3).toInt();
            m_XY = NULL;
            m_Key = 0;
        }
    } else {
        m_XY = NULL;
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

void BSplineVisDialog::showEvent(QShowEvent* event) {

}

void BSplineVisDialog::closeEvent(QCloseEvent* event) {
    
}


