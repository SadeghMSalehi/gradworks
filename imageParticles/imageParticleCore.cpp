//
//  imageParticleCore.cpp
//  imageParticles
//
//  Created by Joohwi Lee on 11/4/12.
//
//

#include "imageParticleTypes.h"
#include "imageParticleCore.h"

class BoundedGradientDescentOptimizer : public itk::RegularStepGradientDescentBaseOptimizer {
public:
    typedef BoundedGradientDescentOptimizer Self;
    typedef RegularStepGradientDescentBaseOptimizer Superclass;
    typedef itk::SmartPointer<Self>       Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    void SetNumberOfParameters(int n) {
        
    }

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(BoundedGradientDescentOptimizer,
                 RegularStepGradientDescentBaseOptimizer);

    /** Advance one step following the gradient direction. */
    virtual void StepAlongGradient(double factor,
                                   const DerivativeType & transformedGradient) {

    }
protected:
    BoundedGradientDescentOptimizer() {}
    virtual ~BoundedGradientDescentOptimizer() {}
};

class LineEntropyCostFunction: public itk::SingleValuedCostFunction {
public:
    typedef LineEntropyCostFunction Self;
    typedef itk::SingleValuedCostFunction Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkTypeMacro(Self, Superclass);
    itkNewMacro(Self);

    typedef double MeasureType;
    typedef Superclass::ParametersType ParametersType;
    typedef Superclass::ParametersValueType ParametersValueType;
    typedef itk::Array<ParametersValueType> DerivativeType;

    void SetNumberOfPoints(int n) {
        m_NumberOfPoints = n;
    }

    void SetNumberOfLambdas(int n) {
        m_NumberOfLambdas = n;
    }

    void SetLowerBounds(int x0) {
        m_LowerBounds = x0;
    }

    void SetUpperBounds(int x1) {
        m_UpperBounds = x1;
    }

    virtual unsigned int GetNumberOfParameters() const {
        return m_NumberOfPoints + m_NumberOfLambdas * 2;
    }

    /**
     * This method returns the value of the cost function corresponding
     * to the specified parameters.    */
    virtual MeasureType GetValue(const ParametersType & parameters) const {
        MeasureType value;
        DerivativeType derivative;
        derivative.SetSize(GetNumberOfParameters());
        GetValueAndDerivative(parameters, value, derivative);
        return value;
    }

    /** This method returns the derivative of the cost function corresponding
     * to the specified parameters.   */
    virtual void GetDerivative(const ParametersType & parameters,
                               DerivativeType & derivative) const {
        MeasureType value;
        GetValueAndDerivative(parameters, value, derivative);
    }


    /** This method returns the value and derivative of the cost function corresponding
     * to the specified parameters    */
    virtual void GetValueAndDerivative(const ParametersType & p,
                                       MeasureType & value,
                                       DerivativeType & derivative) const {
        if (p.GetSize() != m_NumberOfPoints + 2 * m_NumberOfLambdas) {
            return;
        }

        cout << "n Parameters: " << GetNumberOfParameters() << endl;
        MeasureType cost = 0;
        derivative.SetSize(GetNumberOfParameters());
        arma::vec weights;
        weights.zeros(m_NumberOfPoints);

        double si = 3, si2 = si * si;
        for (int i = 0; i < m_NumberOfPoints; i++) {
            for (int j = 0; j < m_NumberOfPoints; j++) {
                if (i == j) {
                    continue;
                }
                double dx = ::abs(p[i] - p[j]);
                double w = exp(-dx * dx / si2);
                if (dx < -3*si || dx > 3*si) {
                    w = 0;
                }
                weights[j] = w;
            }

            double g = arma::sum(weights);

            cost += g;
            derivative[i] = 0;
            if (g <= 0) {
                derivative[i] = 0;
            } else {
                for (int j = 0; j < m_NumberOfPoints; j++) {
                    if (i == j) {
                        continue;
                    }
                    double dx = ::abs(p[i] - p[j]);
                    if (dx > -3*si && dx < 3*si) {
                        derivative[i] -= (p[i] - p[j]) * weights[j] / g;
                    }
                }
                if (p[i] <= m_LowerBounds) {
                    derivative[i] -= p[i + m_NumberOfLambdas];
                } else if (p[i] >= m_UpperBounds) {
                    derivative[i] -= p[i + 2 * m_NumberOfLambdas];
                }
            }
        }

        for (int i = m_NumberOfPoints; i < m_NumberOfLambdas + m_NumberOfPoints; i++) {
            double x = p[i - (m_NumberOfPoints)];
            if (x >= m_LowerBounds) {
                derivative[i] = 0;
                continue;
            }
            cost -= p[i] * (x - m_LowerBounds);
            derivative[i] = - (x - m_LowerBounds);
        }

        for (int i = m_NumberOfPoints + m_NumberOfLambdas; i < m_NumberOfPoints + 2*m_NumberOfLambdas; i++) {
            double x = p[i - (m_NumberOfPoints + m_NumberOfLambdas)];
            if (x <= m_UpperBounds) {
                derivative[i] = 0;
                continue;
            } else {
                cost -= p[i] * (x - m_UpperBounds);
                derivative[i] = - (x - m_UpperBounds);
            }
        }


        value = cost;

        cout << "Parameters: ";
        for (int i = 0; i < m_NumberOfPoints; i++) {
            cout << p[i] << " ";
        }
        cout << endl;
        cout << "Derivatives: " << derivative << endl;

        cout << "Value: " << value << endl;
        //cout << "Derivatives: " << derivative << endl;
    }

protected:
    LineEntropyCostFunction() {
        m_NumberOfPoints = 0;
        m_NumberOfLambdas = 0;
        m_LowerBounds = 0;
        m_UpperBounds = 0;
    }

    virtual ~LineEntropyCostFunction() {
    }
    
private:
    int m_NumberOfPoints;
    int m_NumberOfLambdas;
    int m_LowerBounds;
    int m_UpperBounds;
};

void ImageParticleCore::Run() {
    int N = 5;
    std::srand(100);
    arma::vec points;
    points.randu(N);
    points = points * 100;

    OptimizerType::ParametersType param;
    param.SetSize(3*N);
    for (int i = 0; i < 3*N; i++) {
        if (i < N) {
            param[i] = i * .5;
        } else if (i < 2*N) {
            param[i] = 20;
        } else {
            param[i] = 20;
        }
    }

    LineEntropyCostFunction::Pointer costFunc = LineEntropyCostFunction::New();
    costFunc->SetNumberOfPoints(N);
    costFunc->SetNumberOfLambdas(N);
    costFunc->SetLowerBounds(0);
    costFunc->SetUpperBounds(5);

    OptimizerType::Pointer opti = OptimizerType::New();
    opti->SetCostFunction(costFunc);
    //opti->SetUseUnitLengthGradient(true);
    //opti->SetNumberOfIterations(200);

//    OptimizerType::BoundSelectionType bounds;
//    OptimizerType::BoundValueType lowerBounds, upperBounds;
//
//    bounds.SetSize(3*N);
//    lowerBounds.SetSize(3*N);
//    upperBounds.SetSize(3*N);
//
//    for (int i = 0; i < 3*N; i++) {
//        if (i < N) {
//            bounds[i] = 2;
//        } else {
//            bounds[i] = 0;
//        }
//        lowerBounds[i] = 0;
//        upperBounds[i] = 10;
//    }
//
//    opti->SetBoundSelection(bounds);
//    opti->SetLowerBound(lowerBounds);
//    opti->SetUpperBound(upperBounds);
//    opti->SetTrace(true);

    opti->SetInitialPosition(param);
    opti->StartOptimization();
    cout << opti->GetCurrentPosition() << endl;
}