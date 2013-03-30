//
//  qutils.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/29/13.
//
//

#include "qutils.h"
#include <QFileDialog>
#include <QFileInfo>

QFileManager __fileManager;

void QFileManager::putFile(int key, QString fileName) {
    QFileInfo fileInfo(fileName);
    _sourceDirectories[key] = fileInfo.absolutePath();
}


void QFileManager::putDir(int key, QString dir) {
    _sourceDirectories[key] = dir;
}

QString QFileManager::getDir(int key) {
    return _sourceDirectories[key];
}

QString QFileManager::openFile(int key, QWidget* parent, QString msg, QString filter, QString dir) {
    if (dir.isEmpty()) {
        dir = _sourceDirectories[key];
    }
    std::cout << dir.toStdString() << std::endl;
    QString fileName = QFileDialog::getOpenFileName(parent, msg, dir, filter);
    if (fileName.isEmpty()) {
        return fileName;
    }
    putFile(key, fileName);
    return fileName;
}

QString QFileManager::saveFile(int key, QWidget* parent, QString msg, QString filter, QString dir) {
    if (dir.isEmpty()) {
        dir = _sourceDirectories[key];
    }
    QString fileName = QFileDialog::getSaveFileName(parent, msg, dir, filter);
    if (fileName.isEmpty()) {
        return fileName;
    }
    putFile(key, fileName);
    return fileName;
}
