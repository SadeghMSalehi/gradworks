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
    class ParticleTrace {
    public:
        int GetNumberOfSubject();
        int GetNumberOfTimeSteps(int n);
        int GetMaxId(int n);
        void Clear();
        void Resize(int n);
        void Add(double t, ParticleSubject& subj);
        void Add(double t, ParticleArray& array, int subj = 0);
        void Write(std::ostream& os);
        void Read(std::istream& is);
        static void Read(std::istream& is, ParticleVector& trace);

    private:
        bool AddParticle(Particle& p, int subj = -1);
        
        std::vector<int> m_MaxIds;
        std::vector<int> m_TimeSteps;
        std::vector<Particle> m_BoundingBoxes;
        std::vector<ParticleVector> m_SystemTrace;
    };
}

#endif /* defined(__ParticlesGUI__piParticleTrace__) */
