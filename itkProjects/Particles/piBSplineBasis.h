//
//  piBSplineBasis.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/4/13.
//
//

#ifndef __ParticleGuidedRegistration__piBSplineBasis__
#define __ParticleGuidedRegistration__piBSplineBasis__

#include <iostream>
#include <vector>
#include <cassert>

namespace pi {
    class BSplineBasis {
    public:
        typedef float TReal;
        typedef std::vector<TReal> TRealVector;

    private:
        // knot vectors
        const TReal* _x;
        const int _nx;
        // the index of knot vectors
        const int _i;
        // order (not degree)
        const int _k;

        const TReal _knot;
        const TReal _knotNext;

        BSplineBasis* _lower;
        BSplineBasis* _lowerNext;
        
    public:
        BSplineBasis(TReal* x, int nx, int i, int k)
                :_x(x), _nx(nx), _i(i), _k(k), _knot(x[i]), _knotNext(x[i+k]), _lower(NULL), _lowerNext(NULL) {
            assert(i + k < _nx);
            if (k > 1) {
                _lower = new BSplineBasis(x, nx, i, k - 1);
                if (i + k < _nx) {
                    _lowerNext = new BSplineBasis(x, nx, i + 1, k - 1);
                }
            }
        }
        ~BSplineBasis() {
            if (_lower != NULL) {
                delete _lower;
            }
            if (_lower != NULL) {
                delete _lowerNext;
            }
        }
        inline BSplineBasis& lower() { return *_lower; }
        inline BSplineBasis& lowerNext() { return *_lowerNext; }
        inline TReal lower(double t) { return (*_lower)(t); }
        inline TReal lowerNext(double t) { return (*_lowerNext)(t); }
        inline TReal eval(int i, int k, double t) {
            if (k > _k) {
                return 0;
            }
            if (k == _k) {
                return operator()(t);
            }
            if (i < _i || i > _i + k) {
                return 0;
            }
            if (i == _i) {
                return _lower->eval(i, k, t);
            } else if (_lowerNext != NULL) {
                return _lowerNext->eval(i, k, t);
            }
            return 0;
        }
        
        inline int degree() const { return _k - 1; }
        inline TReal operator()(double t) const {
            if (_k == 1) {
                if (_knot <= t && t < _knotNext) {
                    return 1;
                } else {
                    return 0;
                }
            } else {
                if (t < _x[_i] || t >= _x[_i+_k]) {
                    return 0;
                } else {
                    const double lowerBasis = (*_lower)(t);
                    const double lowerCoeff = (t-_x[_i])/(_x[_i+_k-1]-_x[_i]);
                    if (_lowerNext == NULL) {
                        return lowerCoeff * lowerBasis;
                    }
                    const double lowerNextBasis = (*_lowerNext)(t);
                    const double lowerNextCoeff = (_x[_i+_k]-t)/(_x[_i+_k]-_x[_i+1]);
                    const double value = lowerCoeff * lowerBasis + lowerNextCoeff * lowerNextBasis;
                    return value;
                }
            }
        }
    };
}
#endif /* defined(__ParticleGuidedRegistration__piBSplineBasis__) */