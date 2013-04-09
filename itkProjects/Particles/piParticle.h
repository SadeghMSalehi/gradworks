//
//  piParticle.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/5/13.
//
//

#ifndef ParticleGuidedRegistration_piParticle_h
#define ParticleGuidedRegistration_piParticle_h

#include <iostream>
#include <vector>
#include <vnl/vnl_vector.h>

namespace pi {
    typedef float DataReal;

    // VNL related types
    typedef vnl_vector<DataReal> VNLVector;
    typedef vnl_matrix<DataReal> VNLMatrix;

    class Particle {
    public:
        int subj;
        int idx;
        int label;

        DataReal t;

        // the position x and the transformed point y
        DataReal x[4];
        DataReal y[4];
        DataReal z[4];
        DataReal w[4];


        // the current velocity v and the force f
        DataReal v[4];
        DataReal f[4];

        // the density and pressure of a particle
        DataReal density;
        DataReal pressure;

        bool collisionEvent;

        Particle();
        ~Particle();

        void Zero();
        void Sub(const Particle& p, DataReal* nx);
        void AddForce(DataReal* ff, DataReal alpha = 1);
        DataReal Dist2(const Particle& p);

        Particle& operator=(const Particle& other);
    };
    typedef std::vector<Particle> ParticleVector;

    // utility operator overloading
    std::ostream& operator<<(std::ostream& os, const Particle& par);
    std::istream& operator>>(std::istream& is, Particle& par);
}

#endif
