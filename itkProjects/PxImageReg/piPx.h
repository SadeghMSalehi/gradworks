//
//  piPx.h
//  PxImageReg
//
//  Created by Joohwi Lee on 11/5/13.
//
//

#ifndef __PxImageReg__piPx__
#define __PxImageReg__piPx__

#include <iostream>
#include "piImageDef.h"

namespace pi {
    /// A class represents a particle and implements mathematical operations between particles
    /// @property x position in the subject space
    class Px {
    public:
        double x[__Dim];
        Px() {}
        Px(double v) {
            fordim (k) x[k] = v;
        }
        inline bool isnan() {
            fordim (k) {
                if (x[k] != x[k]) {
                    return true;
                }
            }
            return false;
        }
        /*
         inline Px& operator=(int v) {
         fordim (k) x[k] = v;
         return (*this);
         }
         */
        inline double& operator[](int i) {
            return x[i];
        }
        inline const double& operator[](int i) const {
            return x[i];
        }
        inline Px& operator=(const double v) {
            fordim (k) x[k] = v;
            return (*this);
        }
        inline Px& operator+=(const Px& p) {
            fordim (k) x[k] += p.x[k];
            return (*this);
        }
        inline Px& operator-=(const Px& p) {
            fordim (k) x[k] -= p.x[k];
            return (*this);
        }
        inline double operator*(const Px& p) {
            double d = 0;
            fordim (k) d += x[k] += p.x[k];
            return d;
        }
        inline double dist(const Px& o) {
            return sqrt(dist2(o));
        }
        inline double dist2(const Px& o) {
            double s = 0;
            for (int k = 0; k < __Dim; k++) {
                s += ((x[k] - o.x[k])*(x[k] - o.x[k]));
            }
            return s;
        }
        void dot(const Px& o, const Px& f) {
            fordim (k) x[k] = o[k] * f[k];
        }
        typedef double Elem;
        typedef std::vector<Px> Vector;
    };

    /// A class represents attributes for a particle
    /// @property label gives a membership to a particle
    /// @property bound ???
    class PxA {
    public:
        typedef std::vector<PxA> Vector;

        int label;
        bool bound;

        PxA() { label = 0; }
    };

    // utility operator overloading
    std::ostream& operator<<(std::ostream& os, const Px& par);
    std::ostream& operator<<(std::ostream& os, const Px::Vector& par);
}
#endif /* defined(__PxImageReg__piPx__) */
