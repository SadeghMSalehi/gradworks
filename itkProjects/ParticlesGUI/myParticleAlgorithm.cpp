//
//  myParticleAlgorithm.cpp
//  laplacePDE
//
//  Created by Joohwi Lee on 11/14/12.
//
//

#include "myParticleAlgorithm.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "itkBinaryThresholdImageFilter.h"


namespace surface {


    class OptimizerProgress: public itk::Command {
    private:
        int m_Counter;

    public:
        /** Standard class typedefs. */
        typedef OptimizerProgress Self;
        typedef itk::Command Superclass;
        typedef itk::SmartPointer<Self> Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;

        /** Run-time type information (and related methods). */
        itkTypeMacro(OptimizerProgress, itk::Command)
        ;

        itkNewMacro(OptimizerProgress)
        ;

        void SetAlgorithm(ParticleAlgorithm* algo) {
            m_Algo = algo;
        }

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
            if (++m_Counter % 1 == 0) {
                double value = 0;
                if (dynamic_cast<const LBFGSOptimizerType*>(caller) != NULL) {
                    const LBFGSOptimizerType* opti = dynamic_cast<const LBFGSOptimizerType*>(caller);
                    value = opti->GetCachedValue();
                    //                g_costHistory.push_back(opti->GetCachedValue());
                } else if (dynamic_cast<const FRPROptimizerType*>(caller) != NULL) {
                    const FRPROptimizerType* opti = dynamic_cast<const FRPROptimizerType*>(caller);
                    value = opti->GetValue();
                    //                g_costHistory.push_back(opti->GetValue());
                } else if (dynamic_cast<const GDOptimizerType*>(caller) != NULL) {
                    const GDOptimizerType* opti = dynamic_cast<const GDOptimizerType*>(caller);
                    value = opti->GetValue();
                    //                g_costHistory.push_back(opti->GetValue());
                }
                //            if (realCaller->GetCurrentPosition().GetSize() > 0) {
                //                m_ParametersHistory->push_back(realCaller->GetCurrentPosition());
                //            }
                cout << "Iteration: " << m_Counter << "; Cost: " << value << endl;
                m_Algo->ReportParameters(realCaller->GetCurrentPosition());
            }
        }

    protected:
        OptimizerProgress() {
            m_Counter = 0;
        }
        virtual ~OptimizerProgress() {
        }
    private:
        OptimizerProgress(const Self &);        //purposely not implemented
        void operator=(const Self &); //purposely not implemented
        ParticleAlgorithm* m_Algo;
    };

    class SurfaceEntropyCostFunction: public itk::SingleValuedCostFunction {
    public:
        typedef SurfaceEntropyCostFunction Self;
        typedef itk::SingleValuedCostFunction Superclass;
        typedef itk::SmartPointer<Self> Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;
        typedef itk::Array<float> FloatArrayType;

        itkTypeMacro(Self, Superclass)
        ;itkNewMacro(Self)
        ;

        typedef double MeasureType;
        typedef Superclass::ParametersType ParametersType;
        typedef Superclass::ParametersValueType ParametersValueType;
        typedef itk::Array<ParametersValueType> DerivativeType;

        void SetParticleAlgorithm(ParticleAlgorithm::Pointer algo) {
            m_Algorithm = algo;
        }

        /**
         * Implicit 3d surface defined in label map
         */
        void SetImplicitSurface(LabelType::Pointer label) {
            itk::BinaryThresholdImageFilter<LabelType,LabelType>::Pointer filter = itk::BinaryThresholdImageFilter<LabelType, LabelType>::New();
            filter->SetInput(label);
            filter->SetInsideValue(1);
            filter->SetOutsideValue(0);
            filter->SetLowerThreshold(1);
            filter->SetUpperThreshold(255);
            filter->Update();
            LabelType::Pointer thresholdedLabel = filter->GetOutput();

            DistanceMapFilter::Pointer distanceFilter = DistanceMapFilter::New();
            distanceFilter->SetInput(thresholdedLabel);
            distanceFilter->UseImageSpacingOn();

            std::cout << "Generating distance map ..." << std::flush;
            distanceFilter->Update();
            std::cout << "done" << std::endl;

            DistanceMapFilter::VectorImagePointer vectorImage = distanceFilter->GetVectorDistanceMap();
            ImageType::Pointer distanceImage = distanceFilter->GetOutput();

            m_DistanceImage = distanceImage;
            m_DistanceVector = vectorImage;
        }

        virtual unsigned int GetNumberOfParameters() const {
            return m_Algorithm->GetNumberOfPoints() * m_Algorithm->GetDimensions();
        }

        /** This method returns the value of the cost function corresponding
         * to the specified parameters.    */
        virtual MeasureType GetValue(const ParametersType & parameters) const {
            MeasureType value;
            DerivativeType derivative;
            derivative.SetSize(GetNumberOfParameters());
            GetValueAndDerivative(parameters, value, derivative);
            return value;
        }

        inline ParametersType::ValueType DistanceSquared(const ParametersType& p,  int i, int j) const {
            return (p[i]-p[j])*(p[i]-p[j]) + (p[i+1]-p[j+1])*(p[i+1]-p[j+1]) + (p[i+2]-p[j+2])*(p[i+2]-p[j+2]);
        }

        /** This method returns the derivative of the cost function corresponding
         * to the specified parameters.   */
        virtual void GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const {
            MeasureType value;
            GetValueAndDerivative(parameters, value, derivative);
        }

        virtual void GetValueAndDerivative(const ParametersType & params, MeasureType & value, DerivativeType & derivative) const {
            derivative.SetSize(GetNumberOfParameters());
            derivative.Fill(0);

            const int nDims = m_Algorithm->GetDimensions();
            const int nPoints = m_Algorithm->GetNumberOfPoints();

            FloatArrayType weights;
            weights.SetSize(nPoints + 1);
            value = 0;
            //#pragma omp parallel for
            for (int i = 0; i < nPoints; i++) {
                double pointCost = 0;
                int ii = i  *nDims;
                weights.Fill(0);
                for (int j = 0; j < nPoints; j++) {
                    int jj = j * nDims;
                    if (i == j) {
                        continue;
                    }
                    double d2 = DistanceSquared(params, ii, jj);
                    double weight = 0;
                    if (d2 < m_CutOffDistance * m_CutOffDistance) {
                        weight = exp(-d2 / m_Sigma2);
                    }
                    weights[j] = weight;
                    pointCost += weight;
                }
                double weightSum = weights.sum();
                if (weightSum != 0) {
                    weights /= weightSum;
                }

                // boundary distance computation
                for (int j = 0; j < nPoints; j++) {
                    int jj = j*nDims;
                    for (int k = 0; k < nDims; k++) {
                        derivative[ii + k] -= weights[j] * (params[ii + k] - params[jj + k]);
                    }
                }
                //#pragma omp critical(value)
                {
                    value += pointCost;
                }
            }
        }

    protected:
        SurfaceEntropyCostFunction() {
            m_CutOffDistance = 5;
            m_Sigma2 = 1;

        }

        virtual ~SurfaceEntropyCostFunction() {
        }

    private:
        ImageType::Pointer m_DistanceImage;
        DistanceMapFilter::VectorImagePointer m_DistanceVector;

        SurfaceEntropyCostFunction(const Self&);
        void operator=(const Self&);
        ParticleAlgorithm::Pointer m_Algorithm;
        double m_CutOffDistance;
        double m_Sigma2;
    };


    void ParticleAlgorithm::SetModelScene(vtkPropScene *scene) {
        m_Scene = scene;
        m_Points = dynamic_cast<vtkPointSet*>(m_Scene->FindPolyData("mainSurface"));
        m_Points->Print(cout);
    }


    /**
     * Define where particles are constrained
     */
    void ParticleAlgorithm::SetImageList(ImageContainer::List *imageList) {
        m_ImageList = imageList;
        //    m_ImplicitSurface = image;
        //    DistanceMapFilter::Pointer distanceMapFilter = DistanceMapFilter::New();
        //    distanceMapFilter->SetInput(m_ImplicitSurface);
        //    distanceMapFilter->Update();
        //    m_DistanceValueMap = distanceMapFilter->GetOutput();
        //    m_DistanceVectorMap = distanceMapFilter->GetVectorDistanceMap();
    }



    void ParticleAlgorithm::RunOptimization() {
        if (m_Points == NULL) {
            return;
        }
        
        m_NumberOfPoints = m_Points->GetNumberOfPoints();
        m_InitialParameters.SetSize(m_NumberOfPoints * GetDimensions());
        for (int i = 0; i < m_NumberOfPoints; i++) {
            int j = i * m_Dims;
            double* p = m_Points->GetPoint(i);
            for (int k = 0; k < m_Dims; k++) {
                m_InitialParameters[j + k] = p[k];
            }
        }
        m_ResultPoints = m_Points->NewInstance();
        m_ResultPoints->DeepCopy(m_Points);
        
        
        int nIter = m_Props.GetInt("numberOfIterations", 100);


        GDOptimizerType::Pointer opti = GDOptimizerType::New();
        SurfaceEntropyCostFunction::Pointer costFunc = SurfaceEntropyCostFunction::New();
        costFunc->SetParticleAlgorithm(this);
        costFunc->SetImplicitSurface(m_ImageList->at(0)->GetLabel());
        
        opti->SetCostFunction(costFunc);
        opti->SetInitialPosition(m_InitialParameters);
        opti->SetNumberOfIterations(nIter);
        
        OptimizerProgress::Pointer progress = OptimizerProgress::New();
        progress->SetAlgorithm(this);
        opti->AddObserver(itk::IterationEvent(), progress);
        cout << "Starting optimization ..." << endl;
        opti->StartOptimization();
        cout << "Optimization done ..." << endl;
        OptimizerParametersType result = opti->GetCurrentPosition();
        int j = 0;
        vtkPoints* points = m_ResultPoints->GetPoints();
        for (int i = 0; i < m_NumberOfPoints; i++, j += m_Dims) {
            points->SetPoint(i, result[j], result[j + 1], result[j + 2]);
        }
        m_ResultPoints->SetPoints(points);
    }
    
    vtkPointSet* ParticleAlgorithm::GetResultPoints() const {
        return m_ResultPoints;
    }
    
    void ParticleAlgorithm::ReportParameters(const OptimizerParametersType &params) {
        m_ParameterHistory.push_back(params);
    }
    
    
    
    void ParticleAlgorithm::ContinueOptimization() {
        
    }
    
}