//
//  myParticleKernelFunctions.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 1/10/13.
//
//

#ifndef __ParticlesGUI__myParticleKernelFunctions__
#define __ParticlesGUI__myParticleKernelFunctions__

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include "vnlCommon.h"

namespace my {
    class DefaultKernel {
    public:
        DefaultKernel(double h0) : h(h0), h9(h*h*h*h*h*h*h*h*h) {
        }
        inline double H() const {
            return h;
        }
        double value(double r) const {
            double hr = (h*h-r*r);
            return 315.0/64/M_PI/(h9)*hr*hr*hr;
        }
        void deriv(VNLVector& r, VNLVector& o) const {
            double rr = r.two_norm();
            double coeff = -945.0/32.0/M_PI/h9*(h*h-rr*rr)*(h*h-rr*rr);
            for (int i = 0; i < r.size(); i++) {
                o[i] = r[i] * coeff;
            }
        }
        double laplacian(double rr) const {
            double coeff = -945.0/32.0/M_PI/h9*(h*h-rr*rr)*(h*h-rr*rr);
            return coeff*(h*h-rr*rr)*(3*h*h-7*rr*rr);
        }
    private:
        const double h;
        const double h9;

    };

    class SpikyKernel {
    public:
        SpikyKernel(double h0) : h(h0), pih6(M_PI*h*h*h*h*h*h) {}
        inline double H() const {
            return h;
        }
        double value(double r) const {
            return 15.0/pih6*(h-r)*(h-r)*(h-r);
        }
        // r must not be normalized
        void deriv(VNLVector& r, VNLVector& o) const {
            double rr = r.two_norm();
            double coeff = -45/pih6/rr*(h-rr)*(h-rr);
            for (int i = 0; i < r.size(); i++) {
                o[i] = r[i] * coeff;
            }
        }
        double laplacian(double rr) const {
            return -90.0/pih6/rr*(h-rr)*(h-2*rr);
        }
    private:
        const double h;
        const double pih6;
    };

    class ViscosityKernel {
    public:
        ViscosityKernel(double h0) : h(h0), h2(h*h), h3(h*h*h), twopih3(2*M_PI*h3) {}
        inline double H() const {
            return h;
        }
        double value(double r) const {
            if (r < h && r > 0) {
                return 15/twopih3*(-r*r*r/2/h3 + r*r/h2 + h/2/r - 1);
            } else {
                return 0;
            }
        }
        void deriv(VNLVector& r, VNLVector& o) const {
            double rr = r.two_norm();
            double rr3 = rr*rr*rr;
            double coeff = 15/(twopih3)*(-3*rr/2/h3 + 2/h2 - h/2/rr3);
            for (int i = 0; i < r.size(); i++) {
                o[i] = r[i] * coeff;
            }
        }
        double laplacian(double rr) const {
            return 45.0/M_PI/h3/h3*(h-rr);
        }
    private:
        const double h;
        const double h2;
        const double h3;
        const double twopih3;
    };
}

#endif /* defined(__ParticlesGUI__myParticleKernelFunctions__) */
