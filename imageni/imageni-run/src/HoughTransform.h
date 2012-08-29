/*
 * HoughTransform.h
 *
 *  Created on: Jul 27, 2012
 *      Author: joohwile
 */

#ifndef HOUGHTRANSFORM_H_
#define HOUGHTRANSFORM_H_

#include "Image.h"

namespace cmath {
template<typename S>
class HoughTransform {
public:
    typedef CImage<S> AccumulatorType;
private:
    AccumulatorType* _accmImg;
    int _rRange[2], _aRange[2];
    float _eps;
public:
    HoughTransform(int r0, int r1, int a0, int a1) {
        _rRange[0] = __MIN(r0, r1);
        _rRange[1] = __MAX(r0, r1);
        _aRange[0] = __MIN(a0, a1);
        _aRange[1] = __MAX(a0, a1);
        _accmImg = new AccumulatorType((_rRange[1] - _rRange[0]), (_aRange[1] - _aRange[0]));
        _accmImg->zero();
        _eps = 1e-5;
        cout << "HoughTransform" << endl;
    }
    virtual ~HoughTransform() {
        delete _accmImg;
    }
    void setEpsilon(float eps) {
        _eps = eps;
    }
    void compute(std::vector<int>& x, std::vector<int>& y) {
        cout << "Computing points: " << x.size() << endl;
        for (int s = _rRange[0]; s < _rRange[1]; s++) {
            for (int t = _aRange[0]; t < _aRange[1]; t++) {
                float radT = t * PI / 180.0f;
                for (size_t i = 0; i < x.size(); i++) {
                    int xi = x[i], yi = y[i];
                    if (abs(s - (cos(radT) * xi + sin(radT) * yi)) < _eps) {
                        _accmImg->inc(s, t);
                    }
                }
            }
        }
    }
    AccumulatorType* getAccumulator() {
        return _accmImg;
    }
};

}

#endif /* HOUGHTRANSFORM_H_ */
