/*
 * Gradient.h
 *
 *  Created on: Jul 27, 2012
 *      Author: joohwile
 */

#ifndef GRADIENT_H_
#define GRADIENT_H_

#include "Macros.h"
#include "Image.h"

namespace cmath {

class Gradient {
public:
    template<typename T, typename S>
    static void sobelx(CImage<T>& src, CImage<S>& dst) {
        const int w = src.getWidth();
        const int h = src.getHeight();
        if (w != dst.getWidth() || h != dst.getHeight()) {
            return;
        }
        dst.zero();
        const int stride = src.getStride();
	#pragma omp parallel for
        for (int j = 1; j < h - 1; j++) {
            for (int i = 0; i < w; i++) {
                const T* srcPtr = src.getScanline(j);
                S* dstPtr = dst.getScanline(j);
                dstPtr[i] = (srcPtr[i + stride + 1] + 2 * srcPtr[i + 1] + srcPtr[i - stride + 1]
                                - (srcPtr[i + stride - 1] + 2 * srcPtr[i - 1] + srcPtr[i - stride - 1])) / 2.0f;
#ifdef DEBUG
    			if (isnan(dstPtr[i])) {
    				cout << "NaN detected at dstPtr[i] = " << i << endl;
    			}
#endif
            }
        }
    }

    template<typename T, typename S>
        static void sobely(CImage<T>& src, CImage<S>& dst) {
            const int w = src.getWidth();
            const int h = src.getHeight();
            if (w != dst.getWidth() || h != dst.getHeight()) {
                return;
            }
            dst.zero();
            const int stride = src.getStride();
    #pragma omp parallel for
            for (int j = 1; j < h - 1; j++) {
                for (int i = 0; i < w; i++) {
                    const T* srcPtr = src.getScanline(j);
                    S* dstPtr = dst.getScanline(j);
                    dstPtr[i] = (srcPtr[i + stride + 1] + 2 * srcPtr[i + stride] + srcPtr[i + stride - 1]
                                    - (srcPtr[i - stride - 1] + 2 * srcPtr[i - stride] + srcPtr[i - stride + 1])) / 2.0f;
                }
            }
        }

    template<typename T, typename S>
    static void gradientX(CImage<T>& src, CImage<S>& dst) {
        const int w = src.getWidth();
        const int h = src.getHeight();
        if (w != dst.getWidth() || h != dst.getHeight()) {
            return;
        }
        dst.zero();
#pragma omp parallel for
        for (int j = 1; j < h - 1; j++) {
            for (int i = 0; i < w; i++) {
                const T* srcPtr = src.getScanline(j);
                S* dstPtr = dst.getScanline(j);
                dstPtr[i] = (srcPtr[i + 1] - srcPtr[i - 1]) / 2.0f;
            }
        }
    }

    template<typename T, typename S>
    static void gradientY(CImage<T>& src, CImage<S>& dst) {
        const int w = src.getWidth();
        const int h = src.getHeight();
        if (w != dst.getWidth() || h != dst.getHeight()) {
            return;
        }
        dst.zero();
        const int srcStride = src.getStride();
#pragma omp parallel for
        for (int j = 1; j < h - 1; j++) {
            for (int i = 0; i < w; i++) {
                const T* srcPtr = src.getScanline(j);
                S* dstPtr = dst.getScanline(j);
                dstPtr[i] = (srcPtr[i + srcStride] - srcPtr[i - srcStride]) / 2.0f;
            }
        }
    }

    template<typename T>
    static void gradientMagnitude(CImage<T>& gx, CImage<T>& gy, CImage<T>& gg) {
        const int w = gx.getWidth();
        const int h = gx.getHeight();
        if (w != gy.getWidth() || h != gy.getHeight()) {
            return;
        }
        gg.zero();
#pragma omp parallel for
        forXY(j,h,i,w) {
            gg.set(i, j, sqrt(gx(i, j) * gx(i, j) + gy(i, j) * gy(i, j)));
        }
    }
};

}

#endif /* GRADIENT_H_ */
