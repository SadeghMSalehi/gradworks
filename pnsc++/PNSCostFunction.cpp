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
        double r = parameters[0];

        PNSBase::VectorType deriv(nDims), v(nDims - 1);
        deriv.fill(0);
        for (int j = 1; j < nDims; j++) {
            v[j-1] = parameters[j];
        }
        double tau = m_Tau;

        double cost = 0;
        for (int i = 0; i < nSamples; i++) {
            double dist2 = 0;
            for (int j = 0; j < nDims - 1; j++) {
                double diff = v[j] - m_Data->at(j, i);
                dist2 += (diff * diff);
            }
            double dist = sqrt(dist2);
            cost += ((dist - r) * (dist - r));
            deriv[0] += (2 * (r - dist));
            for (int j = 1; j < nDims; j++) {
                if (dist != 0) {
                    deriv[j] += (2 * (dist - r) * (v[j-1] - m_Data->at(j-1, i)) / dist);
                }
            }
        }
        cost /= nSamples;
        cost +=  tau*(M_PI_2 - r)*(M_PI_2 - r);
        deriv[0] = deriv[0] / nSamples - 2*tau*(M_PI_2 - r);
        for (int j = 0; j < nDims; j++) {
            derivative[j] = deriv[j] / nSamples;
        }

        value = cost;
        return;
    }
}
