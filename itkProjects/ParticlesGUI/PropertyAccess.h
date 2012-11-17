//
//  PropertyAccess.h
//  laplacePDE
//
//  Created by Joohwi Lee on 11/15/12.
//
//

#ifndef __laplacePDE__PropertyAccess__
#define __laplacePDE__PropertyAccess__

#include <iostream>

#include "string"

class QWidget;

class PropertyAccess {
public:
    PropertyAccess() {}
    PropertyAccess(QWidget* widget) : m_widget(widget) { }
    inline void SetWidget(QWidget* widget) { m_widget = widget; }
    double GetDouble(std::string name, double val) const;
    int GetInt(std::string name, int val) const;
    bool GetBool(std::string name, bool val) const;
protected:

private:
    QWidget* m_widget;
};

#endif /* defined(__laplacePDE__PropertyAccess__) */
