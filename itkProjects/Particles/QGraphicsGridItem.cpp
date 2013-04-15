//
//  QGraphicsGridItem.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 12/13/12.
//
//

#include "QGraphicsGridItem.h"
#include "QPainter"
#include "iostream"

using namespace std;

QGraphicsGridItem::QGraphicsGridItem() : QGraphicsItem() {
    m_Pen = QPen(Qt::black, 1);
    m_Pen.setCosmetic(true);
    m_Res = 10;
}

QGraphicsGridItem::~QGraphicsGridItem() {
}

void QGraphicsGridItem::SetPen(QPen pen) {
    m_Pen = pen;
}

void QGraphicsGridItem::SetResolution(int res) {
    m_Res = res;
}

void QGraphicsGridItem::ComputeFromBoundingRect(QRectF rect) {
    VNLMatrix fx, fy;
    fx.set_size(rect.width(), rect.height());
    fy.set_size(rect.width(), rect.height());

    // assume gX and gY are physical coordinate system
    for (int i = 0; i < fx.rows(); i++) {
        for (int j = 0; j < fy.cols(); j++) {
            fx[i][j] = i + rect.x();
            fy[i][j] = j + rect.y();
        }
    }

    SetGrid(fx, fy);

}

void QGraphicsGridItem::SetGrid(VNLMatrix fx, VNLMatrix fy) {
    if (fy.rows() == 0 || fx.rows() == 0) {
        return;
    }
    m_GridX = fx;
    m_GridY = fy;

    m_Bounds.setTop(fy.min_value());
    m_Bounds.setLeft(fx.min_value());
    m_Bounds.setBottom(fy.max_value());
    m_Bounds.setRight(fx.max_value());
    
    prepareGeometryChange();
    update();
}

QRectF QGraphicsGridItem::boundingRect() const {
//    cout << "Bounding Box: " << m_Bounds.left() << ", " << m_Bounds.top() << ", " << m_Bounds.bottom() << ", " << m_Bounds.right() << endl;
    return m_Bounds;
}

void QGraphicsGridItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) {
    const int nRows = m_GridX.rows();
    const int nCols = m_GridX.cols();
    const int nX = (nRows - 1) / m_Res + 1;
    const int nY = (nCols - 1) / m_Res + 1;

    painter->setPen(m_Pen);
    QPointF* px = new QPointF[nX];
    QPointF* py = new QPointF[nY];
    for (int i = 0; i < nRows; i += m_Res) {
        int jj = 0;
        for (int j = 0; j < nCols; j += m_Res, jj++) {
            py[jj].setX(m_GridX[i][j]);
            py[jj].setY(m_GridY[i][j]);
        }
        painter->drawPolyline(py, jj);
        if (jj != nY) {
            std::cout << "Error: ny != jj" << std::endl;
        }
    }
    for (int j = 0; j < nCols; j += m_Res) {
        int ii = 0;
        for (int i = 0; i < nRows; i += m_Res, ii++) {
            px[ii].setX(m_GridX[i][j]);
            px[ii].setY(m_GridY[i][j]);
        }
        painter->drawPolyline(px, ii);
    }
    delete[] px;
    delete[] py;
}