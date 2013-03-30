//
//  qutils.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/29/13.
//
//

#ifndef __ParticleGuidedRegistration__qutils__
#define __ParticleGuidedRegistration__qutils__

#include <iostream>

#include <QWidget>
#include <QHash>
#include <QString>


class QFileManager {
private:
    QHash<int,QString> _sourceDirectories;

public:
    enum FileKind { Single, Source, Target, Label, Image, Fixed, Moving, Transform };

    void putFile(int key, QString fileName);
    void putDir(int key, QString dir);
    QString getDir(int key);
    
    QString openFile(int key, QWidget* parent, QString msg, QString filter = "", QString dir = "");
    QString saveFile(int key, QWidget* parent, QString msg, QString filter = "", QString dir = "");
};

extern QFileManager __fileManager;
extern QHash<QString,QString> __stringHash;

#endif /* defined(__ParticleGuidedRegistration__qutils__) */