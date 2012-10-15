//
//  PNSCostFunction.cpp
//  pnsc++
//
//  Created by Joohwi Lee on 10/12/12.
//
//

#include "PNSCostFunction.h"

/** This method returns the value of the cost function corresponding
 * to the specified parameters.    */
PNSCostFunction::MeasureType PNSCostFunction::GetValue(const ParametersType & parameters) const {
    if (m_Data == NULL) {
        return 0;
    }

    MeasureType value;
    DerivativeType derivative;
    derivative.SetSize(GetNumberOfParameters());
    GetValueAndDerivative(parameters, value, derivative);
    return value;
}

/** This method returns the derivative of the cost function corresponding
 * to the specified parameters.   */
void PNSCostFunction::GetDerivative(const ParametersType& parameters, DerivativeType& derivative) const {
    if (m_Data == NULL) {
        return;
    }

    MeasureType value = 0;
    return GetValueAndDerivative(parameters, value, derivative);
}

void PNSCostFunction::GetValueAndDerivative(const ParametersType & parameters,
                                   MeasureType & value,
                                   DerivativeType & derivative) const {
    if (m_Data == NULL) {
        return;
    }

    if (m_ComputeInEuclideanSpace) {
        int nSamples = m_Data->n_cols;
        int nDims = GetNumberOfParameters();
        double cost = 0;
        double r = parameters[0];

        PNSBase::VectorType deriv(nDims), v(nDims - 1);
        deriv.fill(0);
        for (int j = 1; j < nDims; j++) {
            v[j-1] = parameters[j];
        }
        for (int i = 0; i < nSamples; i++) {
            double dist2 = 0;
            for (int j = 0; j < nDims - 1; j++) {
                double diff = v[j] - m_Data->at(j, i);
                dist2 += (diff * diff);
            }
            double dist = sqrt(dist2);
            cost += (dist - r) * (dist - r);
            deriv[0] += 2 * (r - dist);
            for (int j = 1; j < nDims; j++) {
                if (dist != 0) {
                    deriv[j] += (2 * (dist - r) * (v[j-1] - m_Data->at(j-1, i)) / dist);
                }
            }
        }

        for (int j = 0; j < nDims; j++) {
            derivative[j] = deriv[j];
        }

        cout << "Cost: " << cost << "; Derivative: " << derivative << endl;
        value = cost;
        return;
    }

//    PNSBasic::VectorType normal(parameters[0], parameters[1], parameters[2]);
//    normal = normal.normal();
//    double phi = parameters[3];
//
//    double sum = 0;
//    double derivSum[4] = { 0, 0, 0, 0 };
//    for (int i = 0; i < m_Data->Ncols(); i++) {
//        PNSBasic::VectorType p((*m_Data)[0][i], (*m_Data)[1][i], (*m_Data)[2][i]);
//        double idot = p.dot(normal);
//        if (idot > 1 || idot < -1) {
//            cout << "Range ERROR" << endl;
//        }
//        double diff = acos(idot) - phi;
//        sum += (diff * diff);
//        double idotRes2 = abs(1 - idot*idot);
//        if (idotRes2 < 1e-10) {
//            for (int j = 0; j < 3; j++) {
//                derivSum[j] = 0;
//            }
//        } else {
//            for (int j = 0; j < 3; j++) {
//                derivSum[j] += 2 * diff * ((*m_Data)[i][j] / sqrt(1 - idot*idot));
//            }
//        }
//        derivSum[3] += -2 * diff;
//    }
//
//    for (int i = 0 ; i < derivative.GetSize(); i++) {
//        derivative[i] = derivSum[i];
//    }
//
//    value = sum;
}
