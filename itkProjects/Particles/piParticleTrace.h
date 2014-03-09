//
//  piParticleTrace.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/16/13.
//
//

#ifndef __ParticlesGUI__piParticleTrace__
#define __ParticlesGUI__piParticleTrace__

#include <iostream>
#include "piParticleCore.h"
#include "vector"

// assume particle id ranges from 0 to maxId consecutively
// assume the number of particle per subject are same at every time step for every subject
namespace pi {
    typedef std::vector<ParticleVector> ParticleVectorSeries;

    class ParticleSetSeries {
    public:
        ParticleSetSeries();
        int subjId;
        int maxIdx;
        DataReal lastTime;
        Particle boundingBox;
        ParticleVectorSeries timeSeries;
        bool AppendParticle(int subj, Particle& p, DataReal t);
    };

    class ParticleTrace {
    public:
        ParticleTrace() : isFullTrace(false) {}
        bool isFullTrace;
        void Clear();
        void Resize(int n);
        void Add(DataReal t, ParticleSubject& subj);
        void Add(DataReal t, ParticleArray& array, int subj = 0);
        void Write(std::ostream& os);
        void Read(std::istream& is);
        static void Read(std::istream& is, ParticleVector& trace);
        static void Write(std::ostream& os, ParticleVector& trace);

        // property for whole system
        std::vector<ParticleSetSeries> system;

    private:
        bool AddParticle(Particle& p, int subj = -1);
    };

    ostream& operator<<(ostream& os, ParticleTrace& trace);
}

#endif /* defined(__ParticlesGUI__piParticleTrace__) */
