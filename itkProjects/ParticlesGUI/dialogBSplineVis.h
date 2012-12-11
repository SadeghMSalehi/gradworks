#ifndef __dialogBSplineVis_h__
#define __dialogBSplineVis_h__

#include "ui_dialogBSplineVis.h"
#include "QGraphicsScene"
#include "myImageContainer.h"
#include "vector"

class BSplineVisDialog : public QDialog {
    Q_OBJECT

public:
    BSplineVisDialog(QWidget* parent = NULL);
    virtual ~BSplineVisDialog();

public slots:
    void updateScene();
    
    void on_bspView_mousePressed(QMouseEvent* event);
    void on_bspView_mouseReleased(QMouseEvent* event);
    void on_bspView_mouseMoved(QMouseEvent* event);

    void on_addPairButton_clicked();
    void on_bspViewZoom_sliderMoved(int val);
    void on_updateField_clicked();

protected:
    virtual void showEvent(QShowEvent * event);
    virtual void closeEvent(QCloseEvent * e);
    
    
private:
    SliceType::Pointer m_RefImage;
    DisplacementFieldType::Pointer m_Field;
    Ui::BSplineVisDialog ui;
    QGraphicsScene m_Scene;
    std::vector<QRectF> m_VectorList;
    int m_Key;
    int m_Index;
};
#endif
