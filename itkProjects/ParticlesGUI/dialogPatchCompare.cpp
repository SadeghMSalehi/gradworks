//
// Image Slice compare window
//

#include "dialogPatchCompare.h"
#include "iostream"
#include "mainwindow.h"
#include "myImageContainer.h"

using namespace std;

PatchCompareDialog::PatchCompareDialog(QWidget* parent) {
    ui.setupUi(this);

    ui.leftView->setScene(&m_LeftScene);
    ui.rightView->setScene(&m_RightScene);

    m_Parent = dynamic_cast<MainWindow*>(parent);
}

PatchCompareDialog::~PatchCompareDialog() {

}

void PatchCompareDialog::showEvent(QShowEvent *event) {
    clear();
    
    m_LeftScene.addPixmap(m_Parent->m_ImageList[0]->GetPixmap());
    m_RightScene.addPixmap(m_Parent->m_ImageList[1]->GetPixmap());
}

void PatchCompareDialog::closeEvent(QCloseEvent *e) {
    clear();
}

void PatchCompareDialog::clear() {
    m_LeftScene.clear();
    m_RightScene.clear();
}