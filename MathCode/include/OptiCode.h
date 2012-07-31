/*
 * COpti.h
 *
 *  Created on: Jul 20, 2012
 *      Author: joohwile
 */

#ifndef COPTICODE_H_
#define COPTICODE_H_

#include <string>
#include <iostream>
#include <stdio.h>
#include "Clock.h"
#include "MatrixCode.h"

using namespace std;

namespace MathCode {
template<typename I, typename V, typename G, typename H>
class CFunc {
public:
    typedef I InputType;
    typedef V ValueOutputType;
    typedef G GradientOutputType;
    typedef H HessianOutputType;
    virtual ~CFunc() {
    }
    virtual void position(InputType& x) = 0;
    virtual void value(ValueOutputType& value) = 0;
    virtual void gradient(GradientOutputType& grad) = 0;
    virtual void hessian(HessianOutputType& hessian) = 0;
};

template<typename T>
class CFuncInput {
public:
    virtual ~CFuncInput() {
    }
    virtual T& operator[](int) = 0;
    virtual void negate() = 0;
    virtual void multiply(T);
    virtual void normal(CFuncInput& out) = 0;
    virtual void round(CFuncInput& out) = 0;
    virtual void linearCombination(float, CFuncInput& in, CFuncInput& out);
};

template<typename FuncType>
class CBracketing {
public:
    typedef typename FuncType::InputType InputType;

    FuncType* _fun;
    InputType* _x0;
    InputType* _d;

    void setFunc(FuncType* fun) {
        _fun = fun;
    }
    void setPosition(InputType* x) {
        _x0 = x;
    }
    void setDirection(InputType* d) {
        _d = d;
    }
    float goldenSectionSearch(float a, float b, float c, float tau) {
        InputType bd, xd;
        float x;
        if (c - b > b - a) {
            x = b + 0.38196601125010515179541316563436f * (c - b);
        } else {
            x = b - 0.38196601125010515179541316563436f * (b - a);
        }

        bool term2 = abs(c - a) < tau;
        if (term2) {
            return (c + a) / 2.f;
        }
        _x0->linearCombination(b, *_d, bd);
        _x0->linearCombination(x, *_d, xd);
        _fun->position(xd);
        float fx;
        _fun->value(fx);
        _fun->position(bd);
        float fb;
        _fun->value(fb);
        if (fx < fb) {
            if (c - b > b - a) {
                return goldenSectionSearch(b, x, c, tau);
            } else {
                return goldenSectionSearch(a, x, b, tau);
            }
        } else {
            if (c - b > b - a) {
                return goldenSectionSearch(a, x, b, tau);
            } else {
                return goldenSectionSearch(b, x, c, tau);
            }
        }
    }
};

class CNewtonRhapson {

};

template<typename FuncType>
class CGradientDescent {
public:
    typedef typename FuncType::InputType InputType;
    typedef typename FuncType::GradientOutputType GradientType;

    FuncType* _fun;
    float _eps;
    float _step;
    int _maxIter;

    CGradientDescent(FuncType* fun, float step = 1) :
                    _eps(1e-6), _step(step), _maxIter(1000) {
        _fun = fun;
    }

    inline float abs(float x) {
        return x > 0 ? x : -x;
    }

    int optimize(InputType& x0, InputType& out) {
        CBracketing<FuncType> lineSearch;
        lineSearch.setFunc(_fun);

        CClock timer;
        timer.tick();
        InputType x = x0, xs;
        _fun->position(x0);
        float f;
        _fun->value(f);
        float change = -1;
        for (int i = 0; i < _maxIter; i++) {
            _fun->position(x);
            InputType x1;
            GradientType g;
            _fun->gradient(g);
            g.normal(g);
            bool useBracketing = false;
            if (useBracketing) {
                g.negate();
                lineSearch.setPosition(&x);
                lineSearch.setDirection(&g);
                _step = lineSearch.goldenSectionSearch(0, 2.5f, 5.f, 2);
                x.linearCombination(_step, g, x1);
            } else {
                // rough estimation
                g.round(g);
                x.linearCombination(-_step, g, x1);
            }
            _fun->position(x1);
            float f1;
            _fun->value(f1);
            _DBG_(
                            cout << "iteration: i = " << i << "; x = " << x << " => " << x1 << "; dx = " << g << "; step = " << _step << endl);
            change = f1 - f;
            if (abs(change) < _eps || change > 0) {
                xs = x;
                _DBG_(
                                cout << "solution: xs = " << xs << "; f = " << f << "; f1 = " << f1 << "; df = " << change << endl);
                break;
            } else {
                f = f1;
                x = x1;
            }
        }
        out = xs;
        return timer.tock();
    }
};
}
#endif /* COPTI_H_ */
