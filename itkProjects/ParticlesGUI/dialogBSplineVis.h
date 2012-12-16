#ifndef __dialogBSplineVis_h__
#define __dialogBSplineVis_h__

#include "ui_dialogBSplineVis.h"
#include "QGraphicsScene"
#include "myImageContainer.h"
#include "myImageParticlesAlgorithm.h"
#include "vector"

class ImageParticlesAlgorithm;

class BSplineVisDialog : public QDialog {
    Q_OBJECT

public:
    BSplineVisDialog(QWidget* parent = NULL);
    virtual ~BSplineVisDialog();

    void SetImageParticlesAlgorithm(ImageParticlesAlgorithm::Pointer alg) { m_Algo = alg; }

public slots:
    void updateScene();
    
    void on_bspView_mousePressed(QMouseEvent* event);
    void on_bspView_mouseReleased(QMouseEvent* event);
    void on_bspView_mouseMoved(QMouseEvent* event);

    void on_addPairButton_clicked();
//    void on_bspViewZoom_sliderMoved(int val);
    void on_updateField_clicked();
    void on_copyPointsButton_clicked();

protected:
    virtual void showEvent(QShowEvent * event);
    virtual void closeEvent(QCloseEvent * e);
    
    void CreateGridAndCheckerboards(SliceType::Pointer refImage);
    
private:
//    SliceType::Pointer m_BlackImage;
//    SliceType::Pointer m_WhiteImage;
    SliceType::Pointer m_SrcImage;
    SliceType::Pointer m_DstImage;
    SliceType::Pointer m_WarpedSlice;
    SliceType::Pointer m_WarpedLabelSlice;
    
    SliceType::Pointer m_RefImage;
    SliceType::Pointer m_DetJacobian;
    DisplacementFieldType::Pointer m_Field;
    Ui::BSplineVisDialog ui;
    QGraphicsScene m_Scene;
    std::vector<QRectF> m_VectorList;
    int m_Key;
    int m_Index;

    VNLMatrix gX, gY;
    VNLMatrix tX, tY;

    ImageParticlesAlgorithm::Pointer m_Algo;

};
#endif
