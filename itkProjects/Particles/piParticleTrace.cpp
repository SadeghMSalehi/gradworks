//
//  piParticleTrace.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/16/13.
//
//

#include "piParticleTrace.h"


namespace pi {
    ParticleSetSeries::ParticleSetSeries(): subjId(0), maxIdx(0), lastTime(-1) {
        boundingBox.idx = 0;
    }

    bool ParticleSetSeries::AppendParticle(int subj, Particle& p, DataReal t) {
        bool isNewTime = false;
        if (lastTime != t) {
            // new time series
            timeSeries.resize(timeSeries.size() + 1);
            lastTime = t;
            isNewTime = true;
        }

        ParticleVector& snapshot = timeSeries.back();
        snapshot.push_back(p);
        snapshot.back().t = t;
        snapshot.back().subj = subj;
        maxIdx = ::max(maxIdx, p.idx);
        if (boundingBox.idx == 0) {
            forcopy (p.x, boundingBox.x);
            forcopy (p.x, boundingBox.y);
        } else {
            boundingBox.idx = 1;
            formin (p.x, boundingBox.x, boundingBox.x);
            formax (p.y, boundingBox.y, boundingBox.y);
        }
        return isNewTime;
    }

    void ParticleTrace::Clear() {
        system.clear();
    }
    
    void ParticleTrace::Resize(int n) {
        system.resize(n);
    }


    void ParticleTrace::Add(DataReal t, ParticleSubject& subj) {
        Add(t, subj.m_Particles, subj.m_SubjId);
    }

    void ParticleTrace::Add(DataReal t, ParticleArray& array, int subj) {
        if (subj < 0 || subj > 100) {
            cout << "Skipping possibly too large or negative subject index (" << subj << ")" << endl;
            return;
        }

        if (subj >= system.size()) {
            Resize(subj + 1);
        }

        // assume t is different for each Add() call

        const int n = array.size();
        for (int i = 0; i < n; i++) {
            system[subj].AppendParticle(subj, array[i], t);
        }
    }

    void ParticleTrace::Write(std::ostream& os) {
        const int ns = system.size();
        for (int n = 0; n < ns; n++) {
            for (int t = 0; t < system[n].timeSeries.size(); t++) {
                ParticleVector& snapshot = system[n].timeSeries[t];
                const int np = snapshot.size();
                for (int i = 0; i < np; i++) {
                    Particle& p = snapshot[i];
                    os << p.t << " " << " " << p.subj << " " << p.idx;
                    for4 (k) {
                        os << " " << p.x[k];
                    }
                    os << endl;
                }
            }
        }
        os << endl;
    }

    void ParticleTrace::Read(std::istream& is) {
        system.clear();
        while (is.good()) {
            Particle p;
            is >> p.t;
            if (!is.good()) {
                break;
            }
            is >> p.subj >> p.idx;
            if (p.subj < 0 || p.subj > 100) {
                cout << "subject id is negative or larger than 100" << endl;
                continue;
            }

            if (isFullTrace) {
                is >> p;
            } else {
                for4 (k) {
                    is >> p.x[k];
                }
            }

            if (p.subj >= system.size()) {
                system.resize(p.subj + 1);
            }
            ParticleSetSeries& series = system[p.subj];
            series.AppendParticle(p.subj, p, p.t);
        }
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

    ostream& operator<<(ostream& os, ParticleTrace& trace) {
        cout << "# subjects: " << trace.system.size() << endl;
        if (trace.system.size() == 0) {
            return os;
        }
        cout << "# times: " << trace.system[0].timeSeries.size() << endl;
        if (trace.system[0].timeSeries.size() > 0) {
            cout << "# particles: " << trace.system[0].timeSeries[0].size() << endl;
        }
        return os;
    }
}