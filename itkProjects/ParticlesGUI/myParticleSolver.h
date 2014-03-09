//
//  myParticleSolver.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/29/12.
//
//

#ifndef __ParticlesGUI__myParticleSolver__
#define __ParticlesGUI__myParticleSolver__

#include <iostream>
#include "vnlCommon.h"

namespace my {
    void runODETest();

    class LinearODE {
    public:
        LinearODE();
        LinearODE(const LinearODE&);

        // functor for integration
        void operator()(const VNLVector &x, VNLVector& dydt, const double t);
        // functor for observer
        void operator()(const VNLVector &x, const double t);
        
        void SetInitial(const VNLVector& x);
        void SetMatrix(const VNLMatrix& A);
        void Integrate(double t0, double t1, double dt);
        VNLVector& GetResult();

    private:
        VNLMatrix m_A;
        VNLVector m_x;
    };

}

#endif /* defined(__ParticlesGUI__myParticleSolver__) */