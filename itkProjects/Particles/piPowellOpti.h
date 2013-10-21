//
//  piPowellOpti.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 10/20/13.
//
//

#ifndef __ParticleGuidedRegistration__piPowellOpti__
#define __ParticleGuidedRegistration__piPowellOpti__

#include <iostream>
#include <vector>
#include "newuoa.hh"
#include <nlopt.hpp>

namespace pi {
    typedef std::vector<double> PowellParams;

    template <class F>
    class PowellOpti {
    public:

        double minimizeNEWUOA(F& func, PowellParams& initial, int rb = 1e7);
    };

    template <class F>
    double PowellOpti<F>::minimizeNEWUOA(F& func, PowellParams& initial, int rb) {
        return min_newuoa<double,F>(initial.size(), &initial[0], func, rb);
    }
}
#endif /* defined(__ParticleGuidedRegistration__piPowellOpti__) */
