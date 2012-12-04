#ifndef __dialog_Patch_Compare__
#define __dialog_Patch_Compare__

#include <QDialog>
#include <QCloseEvent>
#include <QShowEvent>
#include <QGraphicsScene>

#include "ui_dialogPatchCompare.h"


class MainWindow;

class PatchCompareDialog: public QDialog {
    Q_OBJECT

public:
    PatchCompareDialog(QWidget* parent);
    ~PatchCompareDialog();

protected:
    virtual void showEvent(QShowEvent * event);
    virtual void closeEvent(QCloseEvent * e);
    
private:
    QGraphicsScene m_LeftScene;
    QGraphicsScene m_RightScene;
    
    Ui::Dialog ui;
    MainWindow* m_Parent;

    void clear();
};

#endif