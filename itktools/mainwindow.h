#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui/QMainWindow>
#include "QGraphicsScene"
#include "QtConcurrentRun"
#include <QFutureWatcher>

#include "ui_mainwindow.h"

#include "itkMyCore.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

public slots:
    void sayHello(void);
    void on_loadTransformButton_clicked();
    void on_saveTransformButton_clicked();    
    void on_runRegistrationButton_clicked();
    void on_actionOpenSource_triggered();
    void on_actionOpenTarget_triggered();
    void on_actionOpenLabel_triggered();
    void on_actionShowSource_triggered(bool c) {
        if (c) {
            ui.actionShowTarget->setChecked(false);
        }
        drawImage();
    }
    void on_actionShowTarget_triggered(bool c) {
        if (c) {
            ui.actionShowSource->setChecked(false);
        }
        drawImage();
    }
    void on_opacityDial_valueChanged() {
        _core.SetLabelOpacity(ui.opacityDial->value());
        drawImage();
    }
    void on_actionShowLabel_triggered() {
        drawImage();
    }
    void on_actionZoomIn_triggered() {
        ui.graphicsView->scale(1.5, 1.5);
    }
    void on_actionZoomOut_triggered() {
        ui.graphicsView->scale(2/3., 2/3.);
    }
    void on_sliceSlider_valueChanged() {
        moveSlice();
    }

    void on_currentRegistrationStep_valueChanged(int value) {
        ui.applyTransformCheck->setChecked(true);
        on_applyTransformCheck_stateChanged(Qt::Checked);
    }
    void on_registrationFinished() {
        if (_core._registrationAlgorithm.IsNull()) {
            cout << "Null registration algorithm ..." << endl;
        }
        ui.toolBox->setCurrentIndex(2);
        int numberOfIterations = _registrationWatcher->result().size();
        on_transformAvailable(numberOfIterations);
        delete _registrationWatcher;
        _registrationWatcher = NULL;
    }
    
    void on_transformAvailable(int numberOfIterations) {
        ui.maxRegistrationSteps->setText(QString("%1").arg(numberOfIterations));
        ui.runRegistrationButton->setEnabled(true);
        ui.currentRegistrationStep->setValue(numberOfIterations);
        ui.currentRegistrationStep->setMaximum(numberOfIterations);
        on_currentRegistrationStep_valueChanged(numberOfIterations);
    }

    void on_applyTransformCheck_stateChanged(int check);

    // when slice is changed
    void on_xySlice_toggled() {
        changeSliceView(2);
    }
    void on_yzSlice_toggled() {
        changeSliceView(0);
    }
    void on_zxSlice_toggled() {
        changeSliceView(1);
    }

    void changeSliceView(int dir);


    void loadDefaults();
    void exit(void);

public:
    void moveSlice();
    void drawImage(bool force = true);



private:
    int _viewingDir;
    Ui::MainWindowClass ui;
    QGraphicsScene _scene;
    itkMyCore _core;
    QFutureWatcher<RegistrationMethod::TransformHistoryType> *_registrationWatcher;

};

#endif // MAINWINDOW_H
