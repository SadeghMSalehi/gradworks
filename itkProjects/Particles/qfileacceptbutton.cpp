//
//  qfileacceptbutton.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/25/13.
//
//

#include "qfileacceptbutton.h"
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QMimeData>
#include <QList>
#include <QUrl>

typedef QList<QUrl> UrlList;

using namespace std;

void QFileAcceptButton::dragEnterEvent(QDragEnterEvent *event) {
    if (event->mimeData()->hasUrls()) {
        UrlList urls = event->mimeData()->urls();
        UrlList::ConstIterator iter;
        for (iter = urls.constBegin(); iter != urls.constEnd(); iter++) {
            const QUrl& url = (*iter);
            if (url.scheme() != "file") {
                event->setAccepted(false);
                return;
            }
        }
        event->setAccepted(true);
        event->acceptProposedAction();
    }
}

void QFileAcceptButton::dropEvent(QDropEvent *event) {
    QList<QString> fileNames;
    if (event->mimeData()->hasUrls()) {
        UrlList urls = event->mimeData()->urls();
        const QUrl& url = urls[0];
        if (url.scheme() != "file") {
            return;
        }
        event->acceptProposedAction();
        QString filePath = url.path();
        emit fileDropped(filePath);
    }
}