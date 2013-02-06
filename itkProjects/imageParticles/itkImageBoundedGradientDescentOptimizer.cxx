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
#ifndef _itkImageBoundedGradientDescentOptimizer_hxx
#define _itkImageBoundedGradientDescentOptimizer_hxx

#include "itkImageBoundedGradientDescentOptimizer.h"
#include "itkBresenhamLine.h"
#include "armadillo"

namespace itk {
    static ImageType::IndexType vec2idx(arma::vec prevPoint) {
        ImageType::IndexType prevIdx;
        for (unsigned int i = 0; i < ImageType::GetImageDimension(); i++) {
            prevIdx[i] = prevPoint[i];
        }
        return prevIdx;
    }

    static arma::vec idx2vec(ImageType::IndexType prevIdx) {
        arma::vec prevPoint;
        prevPoint.zeros(ImageType::GetImageDimension());
        for (unsigned int i = 0; i < ImageType::GetImageDimension(); i++) {
            prevPoint[i] = prevIdx[i];
        }
        return prevPoint;
    }

    arma::vec ImageBoundedGradientDescentOptimizer::ProjectionOntoBoundary(arma::vec pos) {
        arma::vec newpos;
        newpos.zeros(pos.n_elem);
        DistanceVectorImageType::IndexType idx;
        for (unsigned int i = 0; i < ImageType::GetImageDimension(); i++) {
            idx[i] = pos[i];
        }
        DistanceVectorImageType::PixelType offset = m_BoundaryDistanceMap->GetPixel(idx);
        for (unsigned int i = 0; i < ImageType::GetImageDimension(); i++) {
            newpos[i] = idx[i] + offset[i];
        }
        return newpos;
    }

    arma::vec ImageBoundedGradientDescentOptimizer::ComputeTangentialMove(arma::vec pos, arma::vec force, double factor) {
        const static int pointDimension = ImageType::GetImageDimension();
        ImageType::IndexType prevIdx = vec2idx(pos);
        GradientImageType::PixelType boundaryGradient = m_BoundaryGradient->GetPixel(prevIdx);
        arma::vec normal;
        normal.zeros(pointDimension);
        for (int i = 0; i < pointDimension; i++) {
            normal[i] = boundaryGradient[i];
        }
        normal = normal / arma::norm(normal, 2);
        arma::vec tangentialForce = force - arma::dot(normal, force) * normal;
        //tangentialForce.print("Tangential Force = ");

        arma::vec newPos;
        newPos.zeros(pointDimension);
        double magTangentialForce = arma::norm(tangentialForce, 2);
        if (magTangentialForce != 0) {
            tangentialForce = tangentialForce / magTangentialForce;
        }
        
        for (int i = 0; i < pointDimension; i++) {
            newPos[i] = pos[i] + tangentialForce[i] * factor;
        }
        //newPos.print("New Position = ");
        arma::vec constrainedPos = ProjectionOntoBoundary(newPos);
        //constrainedPos.print("Constrained Position = ");
        return constrainedPos;
    }

    arma::vec ImageBoundedGradientDescentOptimizer::SearchClosestBoundary(arma::vec prevPoint) {

    }

    /**
     * Advance one Step following the gradient direction
     * This method will be overrided in non-vector spaces
     */
    void ImageBoundedGradientDescentOptimizer::StepAlongGradient(double factor,
                                                                 const DerivativeType & transformedGradient) {
        itkDebugMacro(
                      << "factor = " << factor << "  transformedGradient= " << transformedGradient);

        const static int pointDimension = ImageType::GetImageDimension();
        const unsigned int spaceDimension = m_CostFunction->GetNumberOfParameters();

        ParametersType newPosition(spaceDimension);
        ParametersType currentPosition = this->GetCurrentPosition();

        for (unsigned int j = 0; j < spaceDimension; j += pointDimension) {
            ImageType::IndexType prevIdx, newIdx, boundaryIdx;
            for (int i = 0; i < pointDimension; i++) {
                newPosition[j+i] = currentPosition[j+i] + transformedGradient[j+i] * factor;
                prevIdx[i] = ::round(currentPosition[j+i]);
                newIdx[i] = ::round(newPosition[j+i]);
            }
            if (m_BoundaryInterpolator->EvaluateAtIndex(prevIdx) == 1 && m_BoundaryInterpolator->EvaluateAtIndex(newIdx) < 1) {
                // the point has reached to the boundary
                // stop the point to the closest boundary point
                itk::BresenhamLine<2> lineBuilder;
                itk::BresenhamLine<2>::IndexArray idxArray = lineBuilder.BuildLine(prevIdx, newIdx);
                for (int i = 0; i < (int) idxArray.size(); i++) {
                    if (m_BoundaryInterpolator->EvaluateAtIndex(idxArray[i]) < 1) {
                        boundaryIdx = idxArray[i];
                        break;
                    }
                }
                for (int i = 0; i < pointDimension; i++) {
                    newPosition[j+i] = boundaryIdx[i];
                }
            } else if (m_BoundaryInterpolator->EvaluateAtIndex(prevIdx) < 1) {
                // if the previous point was on the boundary
                arma::vec force, pos, newpos;
                force.zeros(pointDimension), pos.zeros(pointDimension);
                for (int i = 0; i < pointDimension; i++) {
                    force[i] = transformedGradient[j+i];
                    pos[i] = prevIdx[i];
                }
                newpos = ComputeTangentialMove(pos, force, factor);
                for (int i = 0; i < pointDimension; i++) {
                    newPosition[j+i] = newpos[i];
                }
            }
        }
        
        itkDebugMacro(<< "new position = " << newPosition);
        this->SetCurrentPosition(newPosition);
    }
} // end namespace itk

#endif
