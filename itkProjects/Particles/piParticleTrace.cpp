//
//  piParticleTrace.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/16/13.
//
//

#include "piParticleTrace.h"


namespace pi {

    void ParticleTrace::Add(double t, ParticleArray& array) {
        const int n = array.size();
        for (int i = 0; i < n; i++) {
            m_Trace.push_back(array[i]);
            m_Trace.back().t = t;
            m_Trace.back().idx = i;
        }
    }

    void ParticleTrace::Write(std::ostream& os) {
        const int n = m_Trace.size();
        for (int i = 0; i < n; i++) {
            Particle& p = m_Trace[i];
            os << p.t << " " << p.idx;
            fordim (k) {
                os << " " << p.x[k];
            }
            os << endl;
        }
        os << endl;
    }

    void ParticleTrace::Read(std::istream& is) {
        Read(is, m_Trace);
    }
    
    void ParticleTrace::Read(std::istream& is, ParticleVector& trace) {
        m_Trace.clear();
        while (is.good()) {
            Particle p;
            is >> p.t;
            if (!is.good()) {
                break;
            }
            is >> p.idx;
            fordim (k) {
                is >> p.x[k];
            }
            trace.push_back(p);
        }
    }
}