//
//  piParticleTrace.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/16/13.
//
//

#include "piParticleTrace.h"


namespace pi {
    int ParticleTrace::GetNumberOfSubject() {
        return m_SystemTrace.size();
    }

    int ParticleTrace::GetNumberOfTimeSteps(int n) {
        return m_TimeSteps[n];
    }

    int ParticleTrace::GetMaxId(int n) {
        return m_MaxIds[n];
    }

    void ParticleTrace::Clear() {
        m_MaxIds.clear();
        m_TimeSteps.clear();
        m_BoundingBoxes.clear();
        m_SystemTrace.clear();
    }
    
    void ParticleTrace::Resize(int n) {
        m_MaxIds.resize(n);
        m_TimeSteps.resize(n);
        m_BoundingBoxes.resize(n);
        m_SystemTrace.resize(n);
    }
    
    void ParticleTrace::Add(double t, ParticleSubject& subj) {
        Add(t, subj.m_Particles, subj.m_SubjId);
    }

    void ParticleTrace::Add(double t, ParticleArray& array, int subj) {
        if (subj < 0 || subj > 100) {
            cout << "Skipping possibly too large or negative subject index (" << subj << ")" << endl;
            return;
        }

        if (subj >= GetNumberOfSubject()) {
            Resize(subj + 1);
        }

        // assume t is different for each Add() call
        m_TimeSteps[subj] ++;
        const int n = array.size();
        m_MaxIds[subj] = ::max(n, m_MaxIds[subj]);
        for (int i = 0; i < n; i++) {
            Particle &p = array[i];
            if (AddParticle(p, subj)) {
                m_SystemTrace[subj].back().t = t;
                m_SystemTrace[subj].back().subj = subj;
            }
        }
    }

    void ParticleTrace::Write(std::ostream& os) {
        for (int n = 0; n < m_SystemTrace.size(); n++) {
            const int np = m_SystemTrace[n].size();
            for (int i = 0; i < np; i++) {
                Particle& p = m_SystemTrace[n][i];
                os << p.t << " " << " " << p.subj << " " << p.idx;
                fordim (k) {
                    os << " " << p.x[k];
                }
                os << endl;
            }
        }
        os << endl;
    }

    void ParticleTrace::Read(std::istream& is) {
        m_SystemTrace.clear();
        while (is.good()) {
            Particle p;
            is >> p.t;
            if (!is.good()) {
                break;
            }
            is >> p.subj >> p.idx;
            fordim (k) {
                is >> p.x[k];
            }
            AddParticle(p);
        }
    }

    bool ParticleTrace::AddParticle(Particle& p, int subj) {
        if (subj < 0) {
            subj = p.subj;
        }
        if (subj > 100 || subj < 0) {
            cout << "Skipping possibly too large or negative subject index (" << subj << ")" << endl;
            return false;
        }
        if (m_SystemTrace.size() >= subj) {
            Resize(subj + 1);
        }
        m_SystemTrace[subj].push_back(p);
        m_MaxIds[subj] = ::max(m_MaxIds[subj], p.idx);
        if (m_BoundingBoxes[subj].idx == 0) {
            forcopy (p.x, m_BoundingBoxes[subj].x);
            forcopy (p.x, m_BoundingBoxes[subj].y);
        } else {
            m_BoundingBoxes[subj].idx = 1;
            formin (p.x, m_BoundingBoxes[subj].x, m_BoundingBoxes[subj].x);
            formax (p.y, m_BoundingBoxes[subj].y, m_BoundingBoxes[subj].y);
        }
        return true;
    }

    void ParticleTrace::Read(std::istream& is, ParticleVector& trace) {
        trace.clear();
        while (is.good()) {
            Particle p;
            is >> p.t;
            if (!is.good()) {
                break;
            }
            is >> p.subj >> p.idx;
            fordim (k) {
                is >> p.x[k];
            }
            trace.push_back(p);
        }
    }
}