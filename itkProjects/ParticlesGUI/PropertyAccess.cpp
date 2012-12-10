//
//  PropertyAccess.cpp
//  laplacePDE
//
//  Created by Joohwi Lee on 11/15/12.
//
//

#include "PropertyAccess.h"
#include "QWidget"
#include "QDoubleSpinBox"
#include "QAbstractButton"
#include "QAction"

double PropertyAccess::GetDouble(std::string name, double val) const {
    QDoubleSpinBox* w = m_widget->findChild<QDoubleSpinBox*>(name.c_str());
    if (w == NULL) {
        return val;
    }
    return w->value();
}

int PropertyAccess::GetInt(std::string name, int val) const {
    QSpinBox* w = m_widget->findChild<QSpinBox*>(name.c_str());
    if (w == NULL) {
        return val;
    }
    return w->value();
}

bool PropertyAccess::GetBool(std::string name, bool val) const {
    QAbstractButton* w = m_widget->findChild<QAbstractButton*>(name.c_str());
    if (w == NULL) {
        QAction* a = m_widget->findChild<QAction*>(name.c_str());
        if (a == NULL) {
            return val;
        }
        return a->isChecked();
    }
    return w->isChecked();
}