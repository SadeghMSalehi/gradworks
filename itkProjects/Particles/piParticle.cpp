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


    // constructor
    // set every member variable as zero
    Particle::Particle() {
        Zero();
    }

    Particle::~Particle() {

    }

    void Particle::Zero() {
        t = 0;
        subj = 0;
        idx = 0;
        for4(j) {
            x[j] = y[j] = z[j] = v[j] = f[j] = w[j] = 0;
        }
        density = pressure = 0;
        label = 0;
        collisionEvent = false;
    }

    void Particle::Sub(const Particle& p, DataReal* d) {
        fordim(i) {
            d[i] = x[i] - p.x[i];
        }
    }

    void Particle::AddForce(DataReal* ff, DataReal alpha) {
#if 0
#ifndef BATCH
#if DIMENSIONS == 3
        if (std::abs(ff[0]) > 10 || std::abs(ff[1]) > 10 || std::abs(ff[2]) > 10) {
            cout << "too large force at [" << x[0] << "," << x[1] << "," << x[2] << "]" << endl;
        }
#else
        if (std::abs(ff[0]) > 10 || std::abs(ff[1]) > 10) {
            cout << "too large force: " << flush;
            fordim (k) {
                cout << x[k] << "," << ff[k] << "\t";
            }
            cout << endl;
        }
#endif
#endif
#endif
        fordim(i) {
            f[i] += (alpha * ff[i]);
        }
    }


    DataReal Particle::Dist2(const Particle& p) {
        DataReal d[DIMENSIONS];
        fordim(i) {
            d[i] = x[i] - p.x[i];
        }
        DataReal dist2 = 0;
        fordim(k) {
            dist2 += (d[k]*d[k]);
        }
        return dist2;
    }

    Particle& Particle::operator=(const Particle& other) {
        for4(i) {
            x[i] = other.x[i];
            y[i] = other.y[i];
            v[i] = other.v[i];
            f[i] = other.f[i];
            density = other.density;
            pressure = other.pressure;
        }
        return (*this);
    }

    void createParticles(ParticleVector& particles, int subj, int n) {
        particles.resize(n);
        for (int i = 0; i < particles.size(); i++) {
            particles[i].idx = i;
            particles[i].subj = subj;
        }
    }

}