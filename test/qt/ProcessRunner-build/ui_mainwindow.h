/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created: Wed Oct 9 16:33:31 2013
**      by: Qt User Interface Compiler version 4.8.4
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QFrame>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QMainWindow>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QSplitter>
#include <QtGui/QStatusBar>
#include <QtGui/QTableWidget>
#include <QtGui/QToolBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    QHBoxLayout *horizontalLayout;
    QSplitter *splitter;
    QTableWidget *tableWidget;
    QFrame *frame;
    QPushButton *run;
    QPushButton *stop;
    QWidget *layoutWidget;
    QVBoxLayout *verticalLayout;
    QCheckBox *step1;
    QCheckBox *step2;
    QCheckBox *step3;
    QCheckBox *step4;
    QPushButton *pushButton;
    QPushButton *pushButton_2;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(847, 581);
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        horizontalLayout = new QHBoxLayout(centralWidget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(3, 3, 3, 3);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        splitter = new QSplitter(centralWidget);
        splitter->setObjectName(QString::fromUtf8("splitter"));
        splitter->setOrientation(Qt::Vertical);
        tableWidget = new QTableWidget(splitter);
        if (tableWidget->columnCount() < 3)
            tableWidget->setColumnCount(3);
        QTableWidgetItem *__qtablewidgetitem = new QTableWidgetItem();
        tableWidget->setHorizontalHeaderItem(0, __qtablewidgetitem);
        QTableWidgetItem *__qtablewidgetitem1 = new QTableWidgetItem();
        tableWidget->setHorizontalHeaderItem(1, __qtablewidgetitem1);
        QTableWidgetItem *__qtablewidgetitem2 = new QTableWidgetItem();
        tableWidget->setHorizontalHeaderItem(2, __qtablewidgetitem2);
        if (tableWidget->rowCount() < 4)
            tableWidget->setRowCount(4);
        QTableWidgetItem *__qtablewidgetitem3 = new QTableWidgetItem();
        tableWidget->setVerticalHeaderItem(0, __qtablewidgetitem3);
        QTableWidgetItem *__qtablewidgetitem4 = new QTableWidgetItem();
        tableWidget->setVerticalHeaderItem(1, __qtablewidgetitem4);
        QTableWidgetItem *__qtablewidgetitem5 = new QTableWidgetItem();
        tableWidget->setVerticalHeaderItem(2, __qtablewidgetitem5);
        QTableWidgetItem *__qtablewidgetitem6 = new QTableWidgetItem();
        tableWidget->setVerticalHeaderItem(3, __qtablewidgetitem6);
        tableWidget->setObjectName(QString::fromUtf8("tableWidget"));
        tableWidget->setAlternatingRowColors(true);
        splitter->addWidget(tableWidget);
        frame = new QFrame(splitter);
        frame->setObjectName(QString::fromUtf8("frame"));
        frame->setMinimumSize(QSize(0, 120));
        frame->setMaximumSize(QSize(16777215, 120));
        frame->setFrameShape(QFrame::StyledPanel);
        frame->setFrameShadow(QFrame::Raised);
        run = new QPushButton(frame);
        run->setObjectName(QString::fromUtf8("run"));
        run->setGeometry(QRect(250, 70, 114, 32));
        stop = new QPushButton(frame);
        stop->setObjectName(QString::fromUtf8("stop"));
        stop->setGeometry(QRect(370, 70, 114, 32));
        layoutWidget = new QWidget(frame);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(50, 10, 179, 83));
        verticalLayout = new QVBoxLayout(layoutWidget);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        step1 = new QCheckBox(layoutWidget);
        step1->setObjectName(QString::fromUtf8("step1"));
        step1->setChecked(true);

        verticalLayout->addWidget(step1);

        step2 = new QCheckBox(layoutWidget);
        step2->setObjectName(QString::fromUtf8("step2"));
        step2->setChecked(true);

        verticalLayout->addWidget(step2);

        step3 = new QCheckBox(layoutWidget);
        step3->setObjectName(QString::fromUtf8("step3"));
        step3->setChecked(true);

        verticalLayout->addWidget(step3);

        step4 = new QCheckBox(layoutWidget);
        step4->setObjectName(QString::fromUtf8("step4"));
        step4->setChecked(true);

        verticalLayout->addWidget(step4);

        pushButton = new QPushButton(frame);
        pushButton->setObjectName(QString::fromUtf8("pushButton"));
        pushButton->setGeometry(QRect(250, 10, 114, 32));
        pushButton_2 = new QPushButton(frame);
        pushButton_2->setObjectName(QString::fromUtf8("pushButton_2"));
        pushButton_2->setGeometry(QRect(370, 10, 114, 32));
        splitter->addWidget(frame);

        horizontalLayout->addWidget(splitter);

        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 847, 22));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);

        retranslateUi(MainWindow);
        QObject::connect(run, SIGNAL(clicked()), MainWindow, SLOT(runScript()));
        QObject::connect(stop, SIGNAL(clicked()), MainWindow, SLOT(stopScript()));
        QObject::connect(pushButton, SIGNAL(clicked()), MainWindow, SLOT(addRow()));
        QObject::connect(pushButton_2, SIGNAL(clicked()), MainWindow, SLOT(deleteCurrentRow()));
        QObject::connect(tableWidget, SIGNAL(cellDoubleClicked(int,int)), MainWindow, SLOT(editCell(int,int)));

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MainWindow", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem = tableWidget->horizontalHeaderItem(0);
        ___qtablewidgetitem->setText(QApplication::translate("MainWindow", "#id", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem1 = tableWidget->horizontalHeaderItem(1);
        ___qtablewidgetitem1->setText(QApplication::translate("MainWindow", "#labelmap", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem2 = tableWidget->horizontalHeaderItem(2);
        ___qtablewidgetitem2->setText(QApplication::translate("MainWindow", "#group", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem3 = tableWidget->verticalHeaderItem(0);
        ___qtablewidgetitem3->setText(QApplication::translate("MainWindow", "1", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem4 = tableWidget->verticalHeaderItem(1);
        ___qtablewidgetitem4->setText(QApplication::translate("MainWindow", "2", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem5 = tableWidget->verticalHeaderItem(2);
        ___qtablewidgetitem5->setText(QApplication::translate("MainWindow", "3", 0, QApplication::UnicodeUTF8));
        QTableWidgetItem *___qtablewidgetitem6 = tableWidget->verticalHeaderItem(3);
        ___qtablewidgetitem6->setText(QApplication::translate("MainWindow", "4", 0, QApplication::UnicodeUTF8));
        run->setText(QApplication::translate("MainWindow", "Run", 0, QApplication::UnicodeUTF8));
        stop->setText(QApplication::translate("MainWindow", "Stop", 0, QApplication::UnicodeUTF8));
        step1->setText(QApplication::translate("MainWindow", "Compute Thickness", 0, QApplication::UnicodeUTF8));
        step2->setText(QApplication::translate("MainWindow", "Create SPHARM Mesh", 0, QApplication::UnicodeUTF8));
        step3->setText(QApplication::translate("MainWindow", "Particle Correspondence", 0, QApplication::UnicodeUTF8));
        step4->setText(QApplication::translate("MainWindow", "Perform Analysis", 0, QApplication::UnicodeUTF8));
        pushButton->setText(QApplication::translate("MainWindow", "Add Row", 0, QApplication::UnicodeUTF8));
        pushButton_2->setText(QApplication::translate("MainWindow", "Delete Row", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
