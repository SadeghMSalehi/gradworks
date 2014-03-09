//
//  PNSCostFunction.h
//  pnsc++
//
//  Created by Joohwi Lee on 10/12/12.
//
//

#ifndef __pnsc____PNSCostFunction__
#define __pnsc____PNSCostFunction__

#include <iostream>
#include "itkSingleValuedCostFunction.h"
#include "PNSBase.h"

class PNSCostFunction: public itk::SingleValuedCostFunction {
public:
    typedef PNSCostFunction   Self;
    typedef itk::SingleValuedCostFunction Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(PNSCostFunction, itk::SingleValuedCostFunction);

    itkNewMacro(PNSCostFunction);

    /**  MeasureType typedef.
     *  It defines a type used to return the cost function value. */
    typedef double MeasureType;
    /**  ParametersType typedef.
     *  It defines a position in the optimization search space. */
    typedef Superclass::ParametersType      ParametersType;
    typedef Superclass::ParametersValueType ParametersValueType;

    /** DerivativeType typedef.
     *  It defines a type used to return the cost function derivative.  */
    typedef itk::Array< ParametersValueType > DerivativeType;

    virtual unsigned int GetNumberOfParameters(void) const {
        if (m_Data) {
            if (m_ComputeInEuclideanSpace) {
                return m_Data->n_rows;
            } else {
                return m_Data->n_rows + 1;
            }
        }
        return 0;
    };

    /** This method returns the value of the cost function corresponding
     * to the specified parameters.    */
    virtual MeasureType GetValue(const ParametersType & parameters) const;

    /** This method returns the derivative of the cost function corresponding
     * to the specified parameters.   */
    virtual void GetDerivative(const ParametersType & parameters,
                               DerivativeType & derivative) const;

    /** This method returns the value and derivative of the cost function corresponding
     * to the specified parameters    */
    virtual void GetValueAndDerivative(const ParametersType & parameters,
                                       MeasureType & value,
                                       DerivativeType & derivative) const;
    
    void SetData(PNSBase::MatrixType* data) {
        m_Data = data;
    }

    void SetComputeInEuclideanSpace(bool eSpace) {
        m_ComputeInEuclideanSpace = eSpace;
    }

    void SetTau(double tau) {
    	m_Tau = tau;
    }

protected:
    PNSCostFunction() {
        m_Data = NULL;
        m_ComputeInEuclideanSpace = true;
        m_Tau = 0;
    }
    virtual ~PNSCostFunction() {}

private:
    PNSCostFunction(const Self &); //purposely not implemented
    void operator=(const Self &);           //purposely not implemented

    bool m_ComputeInEuclideanSpace;
    PNSBase::MatrixType* m_Data;
    double m_Tau;
};


#endif /* defined(__pnsc____PNSCostFunction__) */
