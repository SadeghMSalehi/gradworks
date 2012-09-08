//
//  MyFunction.h
//  itkcmds
//
//  Created by Joohwi Lee on 9/7/12.
//
//

#ifndef itkcmds_MyFunction_h
#define itkcmds_MyFunction_h

#include "itkSingleValuedCostFunction.h"


class MyFunction : public itk::SingleValuedCostFunction {
public:
    /** Standard class typedefs. */
    typedef MyFunction   Self;
    typedef itk::SingleValuedCostFunction               Superclass;
    typedef itk::SmartPointer< Self >       Pointer;
    typedef itk::SmartPointer< const Self > ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(MyFunction, itk::SingleValuedCostFunction);

    /** Number of parameters used here */
    virtual unsigned int GetNumberOfParameters(void) const {
        return 2;
    }

    /** This method returns the value of the cost function corresponding
     * to the specified parameters.    */
    virtual MeasureType GetValue(const ParametersType & p) const {
        double x = p[0];
        double y = p[1];
        double s = 1;
        double v = (x*x + y*y - 2) / (s*s) * exp(-(x*x+y*y)/2/(s*s));
        return v;
    };

    /** This method returns the derivative of the cost function corresponding
     * to the specified parameters.   */
    virtual void GetDerivative(const ParametersType & p,
                               DerivativeType & derivative) const {
        derivative.SetSize(2);
        derivative[0] = 2*p[0];
        derivative[1] = 2*p[1];
    }

    /** This method returns the value and derivative of the cost function corresponding
     * to the specified parameters    */
    virtual void GetValueAndDerivative(const ParametersType & parameters,
                                       MeasureType & value,
                                       DerivativeType & derivative) const
    {
        value = this->GetValue(parameters);
        this->GetDerivative(parameters, derivative);
    }

protected:
    MyFunction() {}
    virtual ~MyFunction() {}
private:
    MyFunction(const Self &); //purposely not implemented
    void operator=(const Self &);           //purposely not implemented
    
};


#endif
