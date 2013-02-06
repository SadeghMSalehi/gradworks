#ifndef __ntkSphereMeanComputer_h__
#define __ntkSphereMeanComputer_h__

#include <vtkPoints.h>

#include "itkPowellOptimizer.h"
#include "niralSphericalMeanCostFunction.hpp"

namespace niral {
	
class SphereMeanComputer {
	private:
		int _nPoints;
		vtkPoints* _points;
		double _mean[3];
		bool m_useSphereCoordinate;
		
		void Regulate(double* p) {
			double pa = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
			p[0] /= pa;
			p[1] /= pa;
			p[2] /= pa;
		}
			
	public:
		SphereMeanComputer(int n) : _nPoints(n) {
			_points = vtkPoints::New();
			_points->SetNumberOfPoints(_nPoints);
			_mean[0] = _mean[1] = 1/sqrt(2);
			_mean[2] = 0;
			m_useSphereCoordinate = false;
		};
		
		virtual ~SphereMeanComputer() {
			_points->Delete();
		}
		
		void SphereCoordinateOn() {
			m_useSphereCoordinate = true;
		}
		
		void SphereCoordinateOff() {
			m_useSphereCoordinate = false;
		}
		
		vtkPoints* GetPoints() {
			return _points;
		}
		
		void GetMean(double &x, double &y, double &z) {
			x = _mean[0]; 
			y = _mean[1];
			z = _mean[2];
		}
		
		void ComputeMean() {
			SphericalMeanCostFunction::Pointer costFunc = new SphericalMeanCostFunction;
			costFunc->UnRegister();
			costFunc->SetPoints(_points);
			
			itk::PowellOptimizer::Pointer optimizer = itk::PowellOptimizer::New();
			optimizer->SetCostFunction(costFunc);
			optimizer->SetMaximumIteration(5);
			// optimizer->SetStepTolerance(1e-3);
			
			SphericalMeanCostFunction::ParametersType params(4);
			params[3] = 1;
			
			for (int i = 0; i < 1; i++) {
				params[0] = _mean[0];
				params[1] = _mean[1];
				params[2] = _mean[2];

				optimizer->SetInitialPosition(params);
				optimizer->StartOptimization();
				cout << optimizer->GetStopConditionDescription() << endl;

				_mean[0] = optimizer->GetCurrentPosition()[0];
				_mean[1] = optimizer->GetCurrentPosition()[1];
				_mean[2] = optimizer->GetCurrentPosition()[2];
				
				Regulate(_mean);
			}
			
			cout << setprecision(16) << _mean[0] << ", " << _mean[1] << ", " << _mean[2] << ", " << optimizer->GetCurrentPosition()[3] << endl;
		}
		
		void ComputeMeanInSphereCoordinate() {
			SphericalMeanCostFunction::Pointer costFunc = new SphericalMeanCostFunction;
			costFunc->UnRegister();
			costFunc->SphereCoordinateOn();
			costFunc->SetPoints(_points);
			
			// cout << "Compute mean in sphere coordinate [" << costFunc->GetNumberOfParameters() << "]" << endl;
			
			itk::PowellOptimizer::Pointer optimizer = itk::PowellOptimizer::New();
			optimizer->SetCostFunction(costFunc);
			optimizer->SetMaximumIteration(5);
			// optimizer->SetStepTolerance(1e-3);
			
			SphericalMeanCostFunction::ParametersType params(2);
			
			for (int i = 0; i < 1; i++) {
				double m[3] = { _mean[0], _mean[1], _mean[2] };
				costFunc->GetSphereCoordinate(m, params);
				
				optimizer->SetInitialPosition(params);
				optimizer->StartOptimization();
				
				// cout << optimizer->GetStopConditionDescription() << endl;

				SphericalMeanCostFunction::ParametersType result(3);
				costFunc->GetCartesianCoordinate(optimizer->GetCurrentPosition(), result);
				_mean[0] = result[0];
				_mean[1] = result[1];
				_mean[2] = result[2];
				
				Regulate(_mean);
			}
			
			// cout << setprecision(16) << _mean[0] << ", " << _mean[1] << ", " << _mean[2] << ", " << optimizer->GetCurrentPosition()[3] << endl;
		}
		
		void AddPoint(double x, double y, double z) {
			_points->InsertNextPoint(x,y,z);
		}
		
		void SetPoints(vtkPoints* p) {
			_points = p;
			
			int n = _points->GetNumberOfPoints();
			double sum[3] = { 0 };
			for (int i = 0; i < n; i++) {
				double* pi = _points->GetPoint(i);
				sum[0] += pi[0];
				sum[1] += pi[1];
				sum[2] += pi[2];
			}
			
			_mean[0] = sum[0]/n;
			_mean[1] = sum[1]/n;
			_mean[2] = sum[2]/n;
			
			Regulate(_mean);
		}
};

};

#endif
