/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created: Mon Aug 27 15:38:40 2012
**      by: Qt User Interface Compiler version 4.8.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QHeaderView>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QStatusBar>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_mainwindowClass
{
public:
    QAction *action_New;
    QAction *action_Close;
    QWidget *centralwidget;
    QMenuBar *menubar;
    QMenu *menu_File;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *mainwindowClass)
    {
        if (mainwindowClass->objectName().isEmpty())
            mainwindowClass->setObjectName(QString::fromUtf8("mainwindowClass"));
        mainwindowClass->resize(800, 600);
        action_New = new QAction(mainwindowClass);
        action_New->setObjectName(QString::fromUtf8("action_New"));
        action_Close = new QAction(mainwindowClass);
        action_Close->setObjectName(QString::fromUtf8("action_Close"));
        centralwidget = new QWidget(mainwindowClass);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        mainwindowClass->setCentralWidget(centralwidget);
        menubar = new QMenuBar(mainwindowClass);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 800, 23));
        menu_File = new QMenu(menubar);
        menu_File->setObjectName(QString::fromUtf8("menu_File"));
        mainwindowClass->setMenuBar(menubar);
        statusbar = new QStatusBar(mainwindowClass);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        mainwindowClass->setStatusBar(statusbar);

        menubar->addAction(menu_File->menuAction());
        menu_File->addAction(action_New);
        menu_File->addAction(action_Close);

        retranslateUi(mainwindowClass);

        QMetaObject::connectSlotsByName(mainwindowClass);
    } // setupUi

    void retranslateUi(QMainWindow *mainwindowClass)
    {
        mainwindowClass->setWindowTitle(QApplication::translate("mainwindowClass", "MainWindow", 0, QApplication::UnicodeUTF8));
        action_New->setText(QApplication::translate("mainwindowClass", "&New", 0, QApplication::UnicodeUTF8));
        action_Close->setText(QApplication::translate("mainwindowClass", "&Close", 0, QApplication::UnicodeUTF8));
        menu_File->setTitle(QApplication::translate("mainwindowClass", "&File", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class mainwindowClass: public Ui_mainwindowClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
