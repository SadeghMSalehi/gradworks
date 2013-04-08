//
//  qfileacceptbutton.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/25/13.
//
//

#ifndef __ParticleGuidedRegistration__qfileacceptbutton__
#define __ParticleGuidedRegistration__qfileacceptbutton__

#include <iostream>
#include <QPushButton>

class QFileAcceptButton: public QPushButton {
    Q_OBJECT

signals:
    void fileDropped(QString filePath);

public:
    QFileAcceptButton(QWidget* parent = NULL): QPushButton(parent) {
    }
    
    void dragEnterEvent(QDragEnterEvent* event);
    void dropEvent(QDropEvent* event);

};
#endif /* defined(__ParticleGuidedRegistration__qfileacceptbutton__) */
