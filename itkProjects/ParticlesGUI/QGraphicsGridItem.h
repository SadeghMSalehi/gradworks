//
//  QGraphicsGridItem.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 12/13/12.
//
//

#ifndef __ParticlesGUI__QGraphicsGridItem__
#define __ParticlesGUI__QGraphicsGridItem__

#include <iostream>

#include "QGraphicsItem"
#include "QPen"
#include "vnlCommon.h"

class QGraphicsGridItem : public QGraphicsItem {
public:
    QGraphicsGridItem();
    virtual ~QGraphicsGridItem();

    void SetPen(QPen pen);
    void SetResolution(int res);
    void SetGrid(VNLMatrix fx, VNLMatrix fy);
    
    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);

private:
    QRectF m_Bounds;
    VNLMatrix m_GridX, m_GridY;
    int m_Res;
    QPen m_Pen;
};

#endif /* defined(__ParticlesGUI__QGraphicsGridItem__) */
