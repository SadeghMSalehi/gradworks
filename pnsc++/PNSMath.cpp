//
//  pnsMath.cpp
//  pnsc++
//
//  Created by Joohwi Lee on 10/11/12.
//
//

#include "PNSMath.h"

#include "PNSCostFunction.h"
#include "itkFRPROptimizer.h"
#include "itkSphereOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkCommand.h"
#include "ExpLogTransform.h"


typedef itk::RegularStepGradientDescentOptimizer OptimizerType;

class OptimizerProgress: public itk::Command {
public:
    /** Standard class typedefs. */
    typedef OptimizerProgress Self;
    typedef itk::Command Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(OptimizerProgress, itk::Command);

    itkNewMacro(OptimizerProgress);

    /** Abstract method that defines the action to be taken by the command. */
    virtual void Execute(itk::Object *caller, const itk::EventObject & event) {
        this->Execute((const itk::Object*) caller, event);
    }

    /** Abstract method that defines the action to be taken by the command.
     * This variant is expected to be used when requests comes from a
     * const Object */
    virtual void Execute(const itk::Object *caller, const itk::EventObject & event) {
        const OptimizerType* realCaller = dynamic_cast<const OptimizerType*>(caller);
        if (realCaller == NULL) {
            return;
        }
        cout << realCaller->GetValue() << "; " << realCaller->GetCurrentPosition() << endl;
    }

protected:
    OptimizerProgress() {}
    virtual ~OptimizerProgress() {}
private:
    OptimizerProgress(const Self &);        //purposely not implemented
    void operator=(const Self &); //purposely not implemented
};


void PNSMath::startOptimization(VectorType n0, double phi, double tau) {
    m_Normal = n0;
    m_Phi = phi;

    MatrixType tangentPoints;
    ExpLogTransform expLog;
    expLog.SetTangentPoint(n0);
    expLog.TransformLog(m_Data, tangentPoints);

    tangentPoints.print("Tangent Points: ");

    PNSCostFunction::Pointer costFunc = PNSCostFunction::New();
    costFunc->SetData(&tangentPoints);
    costFunc->SetComputeInEuclideanSpace(true);
    costFunc->SetTau(tau);

    OptimizerType::Pointer opti = OptimizerType::New();
    OptimizerType::ParametersType initialParams;
    initialParams.SetSize(4);
    initialParams[0] = phi;
    initialParams[1] = n0[0];
    initialParams[2] = n0[1];
    initialParams[3] = n0[2];
    OptimizerType::ScalesType scales;
    scales.SetSize(3);
    scales[0] = 1;
    scales[1] = 1;
    scales[2] = 1;
    //scales[3] = 1;
    opti->SetCostFunction(costFunc);
    opti->SetInitialPosition(initialParams);
    //opti->SetUseUnitLengthGradient(true);
    opti->SetScales(scales);
    opti->SetNumberOfIterations(1000);
    opti->SetMaximumStepLength(0.25);

    OptimizerProgress::Pointer progress = OptimizerProgress::New();
    opti->AddObserver(itk::StartEvent(), progress);
    opti->AddObserver(itk::IterationEvent(), progress);
    opti->StartOptimization();


    OptimizerType::ParametersType solution = opti->GetCurrentPosition();

    m_CenterAtTangent.zeros(3);
    m_CenterAtTangent[0] = solution[1];
    m_CenterAtTangent[1] = solution[2];
    m_CenterAtTangent[2] = 1;

    expLog.TransformExp(m_CenterAtTangent, m_Normal);
    m_Phi = solution[0];
}
