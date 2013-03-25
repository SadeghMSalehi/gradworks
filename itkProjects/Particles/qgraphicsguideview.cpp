//
//  QGraphicsGuideView.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/24/13.
//
//

#include "qgraphicsguideview.h"
#include "QMouseEvent"


QGraphicsGuideView::QGraphicsGuideView(QWidget* parent): QGraphicsView(parent) {
    m_selectedItem = NULL;
}

QGraphicsGuideView::~QGraphicsGuideView() {
    QGraphicsView::~QGraphicsView();
}