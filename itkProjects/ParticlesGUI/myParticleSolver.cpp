//
//  myParticleSolver.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/29/12.
//
//

#include "myParticleSolver.h"
#include "boost/numeric/odeint.hpp"

using namespace std;


namespace my {
    LinearODE::LinearODE() {
        m_x.fill(0);
        m_A.fill(0);
    }

    LinearODE::LinearODE(const LinearODE& that) {
        m_x = that.m_x;
        m_A = that.m_A;
    }
    
    void LinearODE::operator()(const VNLVector &x, VNLVector& dydt, const double t) {
        dydt = m_A*x;
    }

    void LinearODE::operator()(const VNLVector &x, const double t) {
        cout << "t: " << t << "; x: " << x << endl;
    }

    void LinearODE::SetInitial(const VNLVector& x) {
        m_x = x;
    }

    void LinearODE::SetMatrix(const VNLMatrix& A) {
        m_A = A;
    }

    VNLVector& LinearODE::GetResult() {
        return m_x;
    }

    void LinearODE::Integrate(double t0, double t1, double dt) {
        boost::numeric::odeint::integrate((*this), m_x , t0 , t1 , 0.1, (*this));
    }

    void runODETest() {
        VNLDiagMatrix A(3);
        A[0] = 1, A[1] = 2, A[2] = 100;
        VNLVector x(3);
        x[0] = x[1] = x[2] = 1;

        LinearODE ode;
        ode.SetInitial(x);
        ode.SetMatrix(A);
        ode.Integrate(0, 1, 0.1);
        cout << "Result: " << ode.GetResult() << endl;
    }
}
