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
    if (m_widget == NULL) return val;
    QDoubleSpinBox* w = m_widget->findChild<QDoubleSpinBox*>(name.c_str());
    if (w == NULL) {
        return val;
    }
    return w->value();
}

int PropertyAccess::GetInt(std::string name, int val) const {
    if (m_widget == NULL) return val;
    QSpinBox* w = m_widget->findChild<QSpinBox*>(name.c_str());
    if (w == NULL) {
        return val;
    }
    return w->value();
}

bool PropertyAccess::GetBool(std::string name, bool val) const {
    if (m_widget == NULL) return val;
    QAbstractButton* w = m_widget->findChild<QAbstractButton*>(name.c_str());
    if (w == NULL) {
        QAction* a = m_widget->findChild<QAction*>(name.c_str());
        if (a == NULL) {
            return val;
        }
        val = a->isChecked();
        // debug for useEnsembleForce
        // action without Checkable flag always return false
//        std::cout << name << ": " << val << std::endl;
        return a->isChecked();
    }
    return w->isChecked();
}