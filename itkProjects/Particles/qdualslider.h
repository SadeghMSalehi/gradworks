//
//  qdualslider.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/21/13.
//
//

#ifndef __ParticleGuidedRegistration__qdualslider__
#define __ParticleGuidedRegistration__qdualslider__

#include <iostream>
#include "QSlider"
#include "QStyle"

class QDualSlider: public QSlider {
    Q_OBJECT

public:
    QDualSlider(QWidget* parent);
    ~QDualSlider();

    int lowValue() { return _low; }
    int highValue() { return _high; }
    void setLowValue(int low) { _low = low; update(); }
    void setHighValue(int high) { _high = high; update(); }

    void setRealMax(double v) { _realMax = v; }
    void setRealMin(double v) { _realMin = v; }
    double realLowValue() { return _low*(_realMax-_realMin)/(maximum()-minimum())+_realMin; }
    double realHighValue() { return _high*(_realMax-_realMin)/(maximum()-minimum())+_realMin; }
    
    virtual void paintEvent(QPaintEvent* event);
    virtual void mousePressEvent(QMouseEvent* event);
    virtual void mouseMoveEvent(QMouseEvent* event);
    virtual void mouseReleaseEvent(QMouseEvent* event);

private:
    int pick(QPoint pt) {
        if (this->orientation() == Qt::Horizontal) {
            return pt.x();
        } else {
            return pt.y();
        }
    }

    int pixelPosToRangeValue(int pos, QStyleOptionSlider& opt);

signals:
    void lowValueChanged(int n);
    void highValueChanged(int n);

private:
    int _low;
    int _high;
    int _clickOffset;
    int _activeSlider;

    double _realMin;
    double _realMax;

    QStyle::SubControl _pressedControl;
    QStyle::SubControl _hoverControl;
};


#endif /* defined(__ParticleGuidedRegistration__qdualslider__) */