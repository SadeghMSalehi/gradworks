//
//  itkOptimizerCommon.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/20/12.
//
//

#ifndef ParticlesGUI_itkOptimizerCommon_h
#define ParticlesGUI_itkOptimizerCommon_h

/**
 * Defines various optimizer types and common parameter type
 *
 */
#include "itkSingleValuedNonLinearOptimizer.h"
#include "itkLBFGSOptimizer.h"
#include "itkFRPROptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

typedef itk::SingleValuedNonLinearOptimizer OptimizerType;
typedef itk::SingleValuedNonLinearOptimizer::ParametersType OptimizerParametersType;
typedef itk::LBFGSOptimizer LBFGSOptimizerType;
typedef itk::FRPROptimizer FRPROptimizerType;
typedef itk::RegularStepGradientDescentOptimizer GDOptimizerType;

#endif
