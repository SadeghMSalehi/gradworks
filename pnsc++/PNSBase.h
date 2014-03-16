//
//  PNSBase.h
//  pnsc++
//
//  Created by Joohwi Lee on 10/12/12.
//
//

#ifndef pnsc___PNSBase_h
#define pnsc___PNSBase_h


#include "MatrixCode.h"
#include <vector>
#include <armadillo>

#include "time.h"
#include "math.h"

class PNSBase {
public:
    typedef arma::mat MatrixType;
    typedef arma::vec VectorType;

    static int GetDimension() {
        return 3;
    }
    
    static void ComputeRotationMatrix(VectorType b, VectorType a, double alpha, MatrixType& rotation) {
        int nDim = b.n_rows;
        b /= arma::norm(b, 2);
        a /= arma::norm(a, 2);

        rotation.eye(nDim, nDim);

        double bTa = arma::dot(b, a);
        if (abs(bTa - 1.) < 1e-15) {
            return;
        }
        if (abs(bTa + 1.) < 1e-15) {
            rotation *= -1;
            return;
        }
        double theta = alpha == 0 ? acos(bTa) : alpha;
        VectorType c(nDim);
        c = b - bTa * a;
        c /= arma::norm(c, 2);
        for (int i = 0; i < nDim; i++) {
            for (int j = 0; j < nDim; j++) {
                rotation.at(j,i) = rotation.at(j,i) + ::sin(theta)*(a[j]*c[i] - c[j]*a[i]) + (::cos(theta) - 1) * (a[j]*a[i] + c[j]*c[i]);
            }
        }
    }

    static void CreateSphereRandoms(int n, double phi, double sigma, VectorType normal, MatrixType& pOut) {
        int nDim = normal.n_elem;

        VectorType randomPhis(n);
        ComputeNormalRandoms(randomPhis, phi, sigma);

        VectorType randomThetas(n);
        ComputeUniformRandoms(randomThetas, 0, 2 * M_PI);

        VectorType northPole(nDim);
        northPole.fill(0);
        northPole[nDim-1] = 1;

        MatrixType rotation;
        ComputeRotationMatrix(northPole, normal, 0, rotation);

        rotation.print("Rotation = ");
        pOut.zeros(nDim, n);
        for (int i = 0; i < n; i++) {
            double theta = randomThetas[i];
            double phi = randomPhis[i];
            VectorType v(nDim), x(nDim);

            v[0] = sin(phi)*cos(theta);
            v[1] = sin(phi)*sin(theta);
            v[2] = cos(phi);
            x = rotation * v;

            pOut.at(0,i) = x[0];
            pOut.at(1,i) = x[1];
            pOut.at(2,i) = x[2];
        }
    }
    static void ComputeNormalRandoms(VectorType &scalars, double mean, double sigma) {
        srand(time(NULL));
        scalars.randn();
        scalars *= sigma;
        scalars += mean;
    }
    static void ComputeUniformRandoms(VectorType &scalars, double m0, double m1) {
        srand(time(NULL));
        scalars.randu();
        scalars *= (m1 - m0);
        scalars += m0;
    }
};

#endif
