//
//  qdualslider.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/21/13.
//
//

#include "qdualslider.h"
#include "QPaintEvent"
#include "QPainter"
#include "QStyleOptionSlider"

using namespace std;

QDualSlider::QDualSlider(QWidget* widget) {
    _pressedControl = QStyle::SC_None;
    _hoverControl = QStyle::SC_None;

    _clickOffset = 0;

    // 0 for the low, 1 for the high, -1 for both
    _activeSlider = 0;
}

QDualSlider::~QDualSlider() {

}

int QDualSlider::pixelPosToRangeValue(int pos, QStyleOptionSlider& opt) {
    QStyle* style = this->style();
    QRect gr = style->subControlRect(QStyle::CC_Slider, &opt, QStyle::SC_SliderGroove);
    QRect sr = style->subControlRect(QStyle::CC_Slider, &opt, QStyle::SC_SliderHandle);

    int sliderMin = gr.x();
    int sliderLength = sr.width();
    int sliderMax = gr.right() - sliderLength + 1;

    if (orientation() != Qt::Horizontal) {
        sliderLength = sr.height();
        sliderMin = gr.y();
        sliderMax = gr.bottom() - sliderLength + 1;
    }
    return style->sliderValueFromPosition(minimum(), maximum(), pos - sliderMin, sliderMax - sliderMin, opt.upsideDown);
}


void QDualSlider::paintEvent(QPaintEvent *event) {
    QPainter painter(this);
    QStyle* style = this->style();

    int values[2] = { _low, _high };
    for (int i = 0; i < 2; i++) {
        QStyleOptionSlider opt;
        this->initStyleOption(&opt);
        if (i == 0) {
            opt.subControls = QStyle::SC_SliderGroove | QStyle::SC_SliderHandle;
        } else {
            opt.subControls = QStyle::SC_SliderHandle;
        }

        if (this->tickPosition() != QSlider::NoTicks) {
            opt.subControls |= QStyle::SC_SliderTickmarks;
        }
        
        if (this->_pressedControl) {
            opt.activeSubControls = this->_pressedControl;
            opt.state |= QStyle::State_Sunken;
        } else {
            opt.activeSubControls = this->_hoverControl;
        }
        
        opt.sliderPosition = values[i];
        opt.sliderValue = values[i];
        style->drawComplexControl(QStyle::CC_Slider, &opt, &painter, this);
    }

}

void QDualSlider::mousePressEvent(QMouseEvent *event) {
    event->accept();

    QStyle* style = this->style();
    Qt::MouseButton button = event->button();

    if (button == Qt::LeftButton) {
        QStyleOptionSlider opt;
        this->initStyleOption(&opt);
        
        this->_activeSlider = -1;
        int values[2] = { _low, _high };
        for (int i = 0; i < 2; i++) {
            opt.sliderPosition = values[i];
            QStyle::SubControl hit = style->hitTestComplexControl(QStyle::CC_Slider, &opt, event->pos());
            if (hit == QStyle::SC_SliderHandle) {
                _activeSlider = i;
                _pressedControl = hit;

                this->triggerAction(QSlider::SliderMove);
                this->setRepeatAction(QSlider::SliderNoAction);
                this->setSliderDown(true);
                break;
            }
        }

        if (_activeSlider < 0) {
            _pressedControl = QStyle::SC_SliderHandle;
            _clickOffset = pixelPosToRangeValue(pick(event->pos()), opt);
            this->triggerAction(QSlider::SliderMove);
            this->setRepeatAction(QSlider::SliderNoAction);
        }
    } else {
        event->ignore();
    }
}

void QDualSlider::mouseMoveEvent(QMouseEvent *event) {
    if (_pressedControl != QStyle::SC_SliderHandle) {
        event->ignore();
        return;
    }

    event->accept();
    QStyleOptionSlider opt;
    this->initStyleOption(&opt);

    int newPos = pixelPosToRangeValue(pick(event->pos()), opt);

    if (_activeSlider < 0) {
        int offset = newPos - _clickOffset;
        _high += offset;
        _low += offset;

        if (_low < minimum()) {
            int diff = minimum() - _low;
            _low += diff;
            _high += diff;
        }
        if (_high > maximum()) {
            int diff = maximum() - _high;
            _low += diff;
            _high += diff;
        }
    } else if (_activeSlider == 0) {
        if (newPos >= _high) {
            newPos = _high - 1;
        }
        _low = newPos;
    } else {
        if (newPos <= _low) {
            newPos = _low + 1;
        }
        _high = newPos;
    }
    _clickOffset = newPos;

    update();

    if (_activeSlider != 1) {
        emit lowValueChanged(newPos);
        emit realLowValueChanged(realLowValue());
    } else if (_activeSlider != 0) {
        emit highValueChanged(newPos);
        emit realHighValueChanged(realHighValue());
    }
}

void QDualSlider::mouseReleaseEvent(QMouseEvent* event) {
    event->accept();
    QStyleOptionSlider opt;
    this->initStyleOption(&opt);

    _pressedControl = _hoverControl;
    _activeSlider = -1;

    emit sliderMoved(lowValue());

    update();
}

 