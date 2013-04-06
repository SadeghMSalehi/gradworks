//
//  piParticle.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/5/13.
//
//

#include "piMacros.h"
#include "piParticle.h"

using namespace std;

namespace pi {
    ostream& operator<<(ostream& os, const Particle& par) {
        for4(k) { os << par.x[k] << " "; }
        for4(k) { os << par.y[k] << " "; }
        for4(k) { os << par.v[k] << " "; }
        for4(k) { os << par.f[k] << " "; }
        os << par.density << " ";
        os << par.pressure << " ";
        return os;
    }

    istream& operator>>(istream& is, Particle& par) {
        for4(k) { is >> par.x[k]; }
        for4(k) { is >> par.y[k]; }
        for4(k) { is >> par.v[k]; }
        for4(k) { is >> par.f[k]; }
        is >> par.density;
        is >> par.pressure;
        return is;
    }

}