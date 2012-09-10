#include "itkSingleValuedCostFunction.h"
#include "vector"

using namespace std;
using namespace itk;

template <typename MetricType>
class MetaMetrics : public SingleValuedCostFunction {
public:
    typedef vector<typename MetricType::Pointer> MetricArrayType;
    typedef vector<int> IntArrayType;

private:
    MetricArrayType _metrics;
    IntArrayType _metricParamterIndex;

public:
    /** Standard class typedefs. */
    typedef MetaMetrics   Self;
    typedef SingleValuedCostFunction               Superclass;
    typedef SmartPointer< Self >       Pointer;
    typedef SmartPointer< const Self > ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(MetaMetrics, SingleValuedCostFunction);

    /**  MeasureType typedef.
     *  It defines a type used to return the cost function value. */
    typedef double MeasureType;

    /**  ParametersType typedef.
     *  It defines a position in the optimization search space. */
    typedef Superclass::ParametersType      ParametersType;
    typedef Superclass::ParametersValueType ParametersValueType;

    /** DerivativeType typedef.
     *  It defines a type used to return the cost function derivative.  */
    typedef Array< ParametersValueType > DerivativeType;

    /** Metric Vector typedef.
     *  It defines a type used to store multiple number of metrics. */
    void AddMetric(typename MetricType::Pointer metric) {
        _metrics.push_back(metric);
        _metricParamterIndex.push_back(metric->GetNumberOfParameters());
    }


    /** Return the number of parameters required to compute
     *  this cost function.
     *  This method MUST be overloaded by derived classes. */
    virtual unsigned int GetNumberOfParameters(void) const  {
        int nParams = 0;
        for (typename MetricArrayType::size_type i = 0 ; i < _metrics.size(); i++) {
            nParams += _metrics[i]->GetNumberOfParameters();
        }
        return nParams;
    }

    /** This method returns the value of the cost function corresponding
     * to the specified parameters.    */
    virtual MeasureType GetValue(const ParametersType & parameters) const {
        MeasureType value = 0;
        int nOffset = 0;
        for (typename MetricArrayType::size_type i = 0; i < _metrics.size(); i++) {
            int nParam = _metrics[i]->GetNumberOfParameters();
            typename MetricType::ParametersType param;
            param.SetSize(nParam);
            for (int j = 0; j < nParam; j++) {
                param[j] = parameters[nOffset + j];
            }
            value += _metrics[i]->GetValue(param);
            nOffset += nParam;
        }
        return value;
    }

    /** This method returns the derivative of the cost function corresponding
     * to the specified parameters.   */
    virtual void GetDerivative(const ParametersType & parameters,
                               DerivativeType & derivative) const {
        derivative.SetSize(this->GetNumberOfParameters());
        int nOffset = 0;
        for (typename MetricArrayType::size_type i = 0; i < _metrics.size(); i++) {
            int nParam = _metrics[i]->GetNumberOfParameters();
            typename MetricType::ParametersType param;
            param.SetSize(nParam);
            for (int j = 0; j < nParam; j++) {
                param[j] = parameters[nOffset + j];
            }
            typename MetricType::DerivativeType deriv;
            _metrics[i]->GetDerivative(param, deriv);
            for (int j = 0; j < nParam; j++) {
                derivative[nOffset + j] = deriv[j];
            }
            nOffset += nParam;
        }
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
    MetaMetrics() {}
    virtual ~MetaMetrics() {}
private:
    MetaMetrics(const Self &); //purposely not implemented
    void operator=(const Self &);           //purposely not implemented
};

