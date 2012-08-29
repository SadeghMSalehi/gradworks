/*
 * HistogramOfGradient.h
 *
 *  Created on: Jul 27, 2012
 *      Author: joohwile
 */

#ifndef HISTOGRAMOFGRADIENT_H_
#define HISTOGRAMOFGRADIENT_H_

#include "Macros.h"
#include "Image.h"

namespace cmath {
template<int N>
class HistogramOfGradient {
private:
    int _width;
    int _height;
    int _n;
public:
    typedef struct bin_s {
        int bin[N];
    } Bin;
    Bin* _d;
    HistogramOfGradient(int w, int h) {
        _width = w;
        _height = h;
        _n = _width * _height;
        _d = new Bin[_n];
        memset(_d, 0, sizeof(Bin) * _n);
    }
    virtual ~HistogramOfGradient() {
        delete[] _d;
    }
    void zero() {
        for (int i = 0; i < _n; i++) {
            memset((void*) &_d[i], 0, sizeof(Bin));
        }
    }
    Bin& operator[](int n) {
        return _d[n];
    }
    void inc(int i, int j, int deg) {
        if (N == 8) {
            int binId = 0;
            if (deg >= 23 && deg < 23 + 45) {
                binId = 1;
            } else if (deg >= 23 + 45 && deg < 23 + 90) {
                binId = 2;
            } else if (deg >= 23 + 90 && deg < 23 + 135) {
                binId = 3;
            } else if (deg >= 23 + 135 && deg < 23 + 180) {
                binId = 4;
            } else if (deg >= 23 + 180 && deg < 23 + 225) {
                binId = 5;
            } else if (deg >= 23 + 225 && deg < 23 + 270) {
                binId = 6;
            } else if (deg >= 23 + 270 && deg < 23 + 315) {
                binId = 7;
            } else if (deg >= 23 + 315 && deg < 23) {
                binId = 0;
            }
            _d[j * _width + i].bin[binId]++;
            cout << "deg: " << deg << "; binId: " << binId << endl;
        } else {
            printf("N more than 8 is not implemented yet!");
        }
    }
    void inc(int i, int j, float y, float x) {
        float deg = atan2(y, x) * 180 / PI;
        inc(i, j, deg);
    }
    void build(CImage<int>& ang) {
        forXY(j, _height, i, _width) {
            inc(i, j, ang(i, j));
        }
    }
    bool createImageRepresentation(CImage<int> &out) {
        if (!out.createEmptyOf(_width * 3, _height * 3)) {
            return false;
        }
        out.zero();
#pragma omp parallel for
        forXY(j,_height, i, _width) {
            const int p = j * _width + i;
            out.set(3 * i + 2, 3 * j + 1, _d[p].bin[0]);
            out.set(3 * i, 3 * j + 1, _d[p].bin[4]);

            out.set(3 * i, 3 * j, _d[p].bin[3]);
            out.set(3 * i + 1, 3 * j, _d[p].bin[2]);
            out.set(3 * i + 2, 3 * j, _d[p].bin[1]);


            out.set(3 * i, 3 * j + 2, _d[p].bin[5]);
            out.set(3 * i + 1, 3 * j + 2, _d[p].bin[6]);
            out.set(3 * i + 2, 3 * j + 2, _d[p].bin[7]);
        }
        return true;
    }
    template<typename T>
    void computeAngle(CImage<T>& gy, CImage<T>& gx, CImage<int>& ang) {
        int w = gy.getWidth();
        int h = gy.getHeight();
        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                T x = gx(i, j);
                T y = gy(i, j);
                int a = round(atan2(y, x) * 180.0f / PI);
                a = (a < 0) ? (a + 360) : a;
                cout << x << ", " << y << "=>" << a << endl;
                ang.set(i, j, a);
            }
        }
    }
    template<typename T>
    static HistogramOfGradient* construct(CImage<T>& gy, CImage<T>& gx) {
        int w = gx.getWidth();
        int h = gx.getHeight();
        HistogramOfGradient<N>* hog = new HistogramOfGradient<N>(w, h);
        CImage<int> ang(w, h);
        hog->computeAngle(gy, gx, ang);
        hog->build(ang);
        return hog;
    }

};
}

#endif /* HISTOGRAMOFGRADIENT_H_ */
