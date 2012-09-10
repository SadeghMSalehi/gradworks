#include "itkCommand.h"
#include "itkRealTimeClock.h"
#include "itkEventObject.h"
#include "iostream"
#include "itkMyFRPROptimizer.h"

class OptimizationReporter : public itk::Command {
public:
    typedef OptimizationReporter Self;
    typedef itk::Command Superclass;
    typedef itk::SmartPointer<Self> Pointer;

    itkTypeMacro(OptimizationReporter, itk::Command);
    itkNewMacro(OptimizationReporter);

    typedef itk::MyFRPROptimizer OptimizerType;

private:

    OptimizationReporter() {
        _updateInterval = 1;
        _numOfIterations = 0;
        _clock = itk::RealTimeClock::New();
    }

    ~OptimizationReporter() {
    }

    itk::RealTimeClock::Pointer _clock;
    itk::RealTimeClock::TimeStampType _lastTime;
    int _numOfIterations;
    int _updateInterval;

public:

    void Execute(itk::Object* caller, const itk::EventObject& event) {
        Execute((const Object*) caller, event);
    }

    void Execute(const itk::Object* object, const itk::EventObject& event) {
        if (object == NULL) {
            std::cout << "Null sender is not processed..." << std::endl;
            return;
        }
        if (typeid(event) == typeid(itk::IterationEvent)) {
            const OptimizerType* opt = dynamic_cast<const OptimizerType*>(object);
            if (++ _numOfIterations % _updateInterval == 0) {
                itk::RealTimeClock::TimeStampType t = _clock->GetTimeInSeconds();
                std::cout << _numOfIterations << "\t" << opt->GetValue() << "\t" << opt->GetCurrentPosition() << "\t" << (t - _lastTime) << " secs" << std::endl;
                _lastTime = t;
            }
        } else if (typeid(event) == typeid(itk::StartEvent)) {
            std::cout << "Optimization has started ..." << std::endl;
            _lastTime = _clock->GetTimeInSeconds();
        } else if (typeid(event) == typeid(itk::EndEvent)) {
            std::cout << "Optimization has finisehd ..." << std::endl;
        }
    }

    void Update() {
        this->Execute((const Object*) NULL, itk::IterationEvent());
    }
};
