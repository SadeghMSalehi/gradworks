//
//  QGraphicsGuideView.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/24/13.
//
//

#ifndef __ParticleGuidedRegistration__QGraphicsGuideView__
#define __ParticleGuidedRegistration__QGraphicsGuideView__

#include <iostream>
#include "QGraphicsView"

class QGraphicsGuideView: public QGraphicsView {
    Q_OBJECT
private:
    QGraphicsItem* m_selectedItem;

public:
    QGraphicsGuideView(QWidget* parent = NULL);
    ~QGraphicsGuideView();

    void mousePressEvent(QMouseEvent* event);
};

#endif /* defined(__ParticleGuidedRegistration__QGraphicsGuideView__) */
