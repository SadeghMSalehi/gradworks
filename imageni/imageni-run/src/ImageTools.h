/*
 * CImageTools.h
 *
 *  Created on: Jul 20, 2012
 *      Author: joohwile
 */

#ifndef CIMAGETOOLS_H_
#define CIMAGETOOLS_H_

#include "MathCode.h"
#include "Image.h"

using namespace cmath;

template<typename T>
class CImageTools {
public:
    typedef Vec4 CRect;

    static CRect computeBoundingBox(Vec3* pts, int n) {
        CRect rect;
        rect[0] = pts[0].x();
        rect[1] = pts[0].y();
        rect[2] = pts[0].x();
        rect[3] = pts[0].y();

        for (int i = 1; i < n; i++) {
            if (pts[i].x() < rect[0]) {
                rect[0] = pts[i].x();
            }
            if (pts[i].y() < rect[1]) {
                rect[1] = pts[i].y();
            }
            if (pts[i].x() > rect[2]) {
                rect[2] = pts[i].x();
            }
            if (pts[i].y() > rect[3]) {
                rect[3] = pts[i].y();
            }
        }

        // set to width and height
        rect[2] = rect[2] - rect[0];
        rect[3] = rect[3] - rect[1];
        return rect;
    }

    /*
    static CImage<T>* rotation(CImage<T>& src, double ang) {
        Mat3 src2dst;
        Vec3 bbox[4];
        //src2dst.rotation(-ang, src.getWidth() / 2, src.getHeight() / 2);
        cout << src2dst << endl;

        //mult(0, 0, 0, bbox[0]);
        //mult(0, src.getWidth(), 0, bbox[1]);
        //mult(0, src.getHeight(), 0, bbox[2]);
        //mult(src.getWidth(), src.getHeight(), 0, bbox[3]);

        for (int i = 0; i < 4; i++) {
            cout << bbox[i] << endl;
        }

        CRect rect = computeBoundingBox(bbox, 4);
        cout << "Bounding box: " << rect << endl;

        CImage<T>* dst = new CImage<T>(rect[2], rect[3]);
        Mat3 dst2src;
        //dst2src.rotation(-ang, src.getWidth() / 2, src.getHeight() / 2);
        for (int j = rect[1]; j < rect[1] + rect[3]; j++) {
            for (int i = rect[0]; i < rect[0] + rect[2]; i++) {
                Vec3 srcPt;
                if (srcPt[0] < 0 || srcPt[0] >= src._width || srcPt[1] < 0 || srcPt[1] >= src._height) {
                    dst2src.mult(i, j, 1, srcPt);
                    (*dst)(i - rect[0], j - rect[1]) = src(srcPt[0], srcPt[1]);
                }
            }
        }
        return dst;
    }
    */

};

#endif /* CIMAGETOOLS_H_ */
