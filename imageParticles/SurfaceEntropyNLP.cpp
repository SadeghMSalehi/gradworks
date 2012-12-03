// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: hs071_nlp.cpp 1324 2008-09-16 14:19:26Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "IpIpoptApplication.hpp"
#include "SurfaceEntropyNLP.h"

// for printf
#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

using namespace Ipopt;

// constructor
SurfaceEntropyNLP::SurfaceEntropyNLP()
{
}

//destructor
SurfaceEntropyNLP::~SurfaceEntropyNLP()
{}

// returns the size of the problem
bool SurfaceEntropyNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                     Index& nnz_h_lag, IndexStyleEnum& index_style)
{
    n = m_NumberOfPoints * m_Dimensions;
    m = m_NumberOfPoints;
    nnz_jac_g = m_NumberOfPoints * m_Dimensions;
    // won't use the hessian of lagrangian
    nnz_h_lag = 0;
    index_style = TNLP::C_STYLE;
    return true;
}

// returns the variable bounds
bool SurfaceEntropyNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                        Index m, Number* g_l, Number* g_u)
{
    // the variables have lower bounds of 1
    for (Index i = 0; i < m_NumberOfPoints * m_Dimensions; i += m_Dimensions) {
        x_l[i] = 0;
        x_l[i+1] = 0;
    }

    // the variables have upper bounds of 5
    for (Index i = 0; i < m_NumberOfPoints * m_Dimensions; i += m_Dimensions) {
        x_u[i] = m_Region[0];
        x_u[i+1] = m_Region[1];
    }

    for (Index i = 0; i < m_NumberOfPoints; i++) {
        g_l[i] = 0;
        g_u[i] = m_Radius[0] * m_Radius[1];
    }
    return true;
}

// returns the initial point for the problem
bool SurfaceEntropyNLP::get_starting_point(Index n, bool init_x, Number* x,
                                           bool init_z, Number* z_L, Number* z_U,
                                           Index m, bool init_lambda,
                                           Number* lambda)
{
    // Here, we assume we only have starting values for x, if you code
    // your own NLP, you can provide starting values for the dual variables
    // if you wish
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);

    for (int i = 0; i < m_NumberOfPoints * m_Dimensions; i++) {
        x[i] = m_InitialPoints(0, i);
    }
    return true;
}

double _si = 10;
static inline double G(double d2, double si) {
    if (d2 > 9*si*si || d2 < -9*si*si) {
        return 0;
    }
    return exp(-(d2) / (2*si*si)) / (sqrt(2*M_PI)*si);
}

// returns the value of the objective function
bool SurfaceEntropyNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
    double gsum = 0;
    for (int i = 0, j = 0; i < m_NumberOfPoints; i++, j += m_Dimensions) {
        for (int k = i, l = j; k < m_NumberOfPoints; k++, l += m_Dimensions) {
            if (i == k) {
                continue;
            }
            double dist2 = (x[j]-x[l])*(x[j]-x[l]) + (x[j+1]-x[l+1])*(x[j+1]-x[l+1]);
            gsum += G(dist2, _si);
        }
    }
    obj_value = gsum * 2;
    return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool SurfaceEntropyNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
    for (int i = 0, l = 0; i < m_NumberOfPoints; i++, l += m_Dimensions) {
        double gx[m_Dimensions] = { 0 };
        double gsum = 0;
        for (int j = 0, m = 0; j < m_NumberOfPoints; j++, m += m_Dimensions) {
            if (i == j) {
                continue;
            }
            double dist2 = (x[l]-x[m])*(x[l]-x[m]) + (x[l+1]-x[m+1])*(x[l+1]-x[m+1]);
            gsum += G(dist2, _si);
        }
        for (int j = 0, m = 0; j < m_NumberOfPoints; j++, m += m_Dimensions) {
            if (i == j) {
                continue;
            }
            double dist2 = (x[l]-x[m])*(x[l]-x[m]) + (x[l+1]-x[m+1])*(x[l+1]-x[m+1]);
            double g = G(dist2, _si);
            for (int k = 0; k < m_Dimensions; k++) {
                if (gsum != 0) {
                    gx[k] += g * (x[l+k] - x[m+k]) / gsum;
                }
            }
        }
        for (int k = 0; k < m_Dimensions; k++) {
            grad_f[l+k] = gx[k] / (_si*_si);
        }
    }
    return true;
}

// return the value of the constraints: g(x)
bool SurfaceEntropyNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    for (int i = 0, j = 0; i < m_NumberOfPoints; i++, j += m_Dimensions) {
        g[i] = (x[j] - m_ImageCenter[0]) * (x[j] - m_ImageCenter[0])
            + (x[j+1] - m_ImageCenter[1]) * (x[j+1] - m_ImageCenter[1]);
    }

    return true;
}

// return the structure or values of the jacobian
bool SurfaceEntropyNLP::eval_jac_g(Index n, const Number* x, bool new_x,
                                   Index m, Index nele_jac, Index* iRow, Index *jCol,
                                   Number* values)
{
    if (values == NULL) {
        int k = 0;
        for (int i = 0; i < m_NumberOfPoints; i ++) {
            for (int j = 0; j < m_Dimensions; j++) {
                iRow[k] = i;
                jCol[k] = j;
                k++;
            }
        }
    }
    else {
        for (int i = 0, j = 0; i < m_NumberOfPoints; i++, j += m_Dimensions) {
            values[j] = 2 * (x[j] - m_ImageCenter[0]);
            values[j+1] = 2 * (x[j+1] - m_ImageCenter[1]);
        }
    }

    return true;
}

//return the structure or values of the hessian
bool SurfaceEntropyNLP::eval_h(Index n, const Number* x, bool new_x,
                               Number obj_factor, Index m, const Number* lambda,
                               bool new_lambda, Index nele_hess, Index* iRow,
                               Index* jCol, Number* values)
{
    return false;
}

void SurfaceEntropyNLP::finalize_solution(SolverReturn status,
                                          Index n, const Number* x, const Number* z_L, const Number* z_U,
                                          Index m, const Number* g, const Number* lambda,
                                          Number obj_value,
                                          const IpoptData* ip_data,
                                          IpoptCalculatedQuantities* ip_cq)
{
    // here is where we would store the solution to variables, or write to a file, etc
    // so we could use the solution.
    m_ResultPoints.zeros(n);
    m_FinalCost = obj_value;

    // For this example, we write the solution to the console
    printf("\n\nSolution of the primal variables, x\n");
    for (Index i=0; i<n; i++) {
        m_ResultPoints[i] = x[i];
        printf("x[%d] = %e\n", i, x[i]);
    }
    
    printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
    for (Index i=0; i<n; i++) {
        printf("z_L[%d] = %e\n", i, z_L[i]);
    }
    for (Index i=0; i<n; i++) {
        printf("z_U[%d] = %e\n", i, z_U[i]);
    }


    printf("\n\nObjective value\n");
    printf("f(x*) = %e\n", obj_value);
    
    printf("\nFinal value of the constraints:\n");
    for (Index i=0; i<m ;i++) {
        printf("g(%d) = %e\n", i, g[i]);
    }
}

bool SurfaceEntropyNLP::StartOptimization() {
    return true;
}