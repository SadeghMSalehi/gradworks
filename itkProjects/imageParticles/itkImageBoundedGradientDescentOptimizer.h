/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkImageBoundedGradientDescentOptimizer_h
#define __itkImageBoundedGradientDescentOptimizer_h

#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkVersor.h"
#include "imageParticleTypes.h"

namespace itk
{
    class ImageBoundedGradientDescentOptimizer: public RegularStepGradientDescentBaseOptimizer {
    public:
        /** Standard class typedefs. */
        typedef ImageBoundedGradientDescentOptimizer         Self;
        typedef RegularStepGradientDescentBaseOptimizer Superclass;
        typedef SmartPointer< Self >                    Pointer;
        typedef SmartPointer< const Self >              ConstPointer;
        typedef itk::Vector<float,2>	VectorType;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(ImageBoundedGradientDescentOptimizer,
                     RegularStepGradientDescentBaseOptimizer);

        void SetBoundaryImage(ImageType::Pointer boundaryImage) {
            m_Boundary = boundaryImage;
            GradientImageFilter::Pointer gradFilter = GradientImageFilter::New();
            gradFilter->SetInput(m_Boundary);
            gradFilter->SetSigma(1);
            gradFilter->Update();
            m_BoundaryGradient = gradFilter->GetOutput();
            m_BoundaryInterpolator = NearestNeighborInterpolatorType::New();
            m_BoundaryInterpolator->SetInputImage(m_Boundary);
            DistanceMapFilter::Pointer distFilter = DistanceMapFilter::New();
            distFilter->SetInput(m_Boundary);
            distFilter->Update();
            m_BoundaryDistanceMap = distFilter->GetVectorDistanceMap();
        }

        arma::vec ProjectionOntoBoundary(arma::vec pos);
        arma::vec ComputeTangentialMove(arma::vec pos, arma::vec force, double factor = 0);
        arma::vec SearchClosestBoundary(arma::vec prevPoint);
        
        /** Advance one step following the gradient direction. */
        virtual void StepAlongGradient(double factor,
                                       const DerivativeType & transformedGradient);



    protected:
        ImageBoundedGradientDescentOptimizer() {}
        virtual ~ImageBoundedGradientDescentOptimizer() {}

    private:
        ImageBoundedGradientDescentOptimizer(const Self &); //purposely not implemented
        void operator=(const Self &);                  //purposely not implemented
        double m_Radius;

        VectorType m_Center;
        ImageType::Pointer m_Boundary;
        GradientImageType::Pointer m_BoundaryGradient;
        DistanceVectorImageType::Pointer m_BoundaryDistanceMap;
        NearestNeighborInterpolatorType::Pointer m_BoundaryInterpolator;

    };
} // end namespace itk

#endif