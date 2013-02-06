#ifndef __niralSphericalMeanCostFunction_hpp__
#define __niralSphericalMeanCostFunction_hpp__

#include "itkSingleValuedCostFunction.h"
#include "vtkPoints.h"

namespace niral {
class SphericalMeanCostFunction : public itk::SingleValuedCostFunction {
	public:
		typedef itk::SingleValuedCostFunction Superclass;
		typedef SphericalMeanCostFunction Self;
		typedef itk::SmartPointer<Self> Pointer;
		typedef itk::SmartPointer<const Self> ConstPointer;
		
		itkFactorylessNewMacro(Self);
		itkTypeMacro(SphericalMeanCostFunction, itk::SingleValuedCostFunction);
		
	private:
		vtkPoints* _points;
		bool m_SphereCoordinate;
		
		SphericalMeanCostFunction(const Self&);
		void operator=(const Self&);
		
	public:
		SphericalMeanCostFunction() : m_SphereCoordinate(false) { };
		virtual ~SphericalMeanCostFunction() {};
		
	public:
		void SetPoints(vtkPoints* p) { 
			_points = p; 
		}
		
		void SphereCoordinateOn() {
			m_SphereCoordinate = true;
		}
		
		void SphereCoordinateOff() {
			m_SphereCoordinate = false;
		}
		
		void GetSphereCoordinate(const ParametersType &c, ParametersType &s) {
			s[0] = atan2(c[1], c[0]);
			s[1] = acos(c[2]);
		}
		
		void GetSphereCoordinate(const double *c, ParametersType &s) {
			s[0] = atan2(c[1], c[0]);
			s[1] = acos(c[2]);
		}		
		
		void GetCartesianCoordinate(const ParametersType &s, ParametersType &c) {
			c[0] = sin(s[1]) * cos(s[0]);
			c[1] = sin(s[1]) * sin(s[0]);
			c[2] = cos(s[1]);
		}
		
		void GetCartesianCoordinate(const double *s, ParametersType &c) {
			c[0] = sin(s[1]) * cos(s[0]);
			c[1] = sin(s[1]) * sin(s[0]);
			c[2] = cos(s[1]);
		}		
		
		unsigned int GetNumberOfParameters() const { return m_SphereCoordinate ? 2 : 4; }
		
		void GetDerivative(const ParametersType &params, DerivativeType &derivs) const {
			return;
		}
		
		MeasureType GetValue(const ParametersType &params) const {
				if (m_SphereCoordinate) {
					return GetValueWithSphereCoordinate(params);
				} else {
					return GetValueWithCartesianCoordinate(params);
				}
		};
		
		double ComputeSquaredAngleSum(const ParametersType &p, double pm) const {
			double angDist2 = 0;

			for (int i = 0; i < _points->GetNumberOfPoints(); i++) {
				double *vi = _points->GetPoint(i);
				double inp = (p[0]*vi[0] + p[1]*vi[1] + p[2]*vi[2]);
				double ang = acos(inp);
				angDist2 += (ang*ang);
			}
			return angDist2;
		};
		
		MeasureType GetValueWithSphereCoordinate(const ParametersType &params) const {
			ParametersType p(3);
			
			p[0] = sin(params[1]) * cos(params[0]);
			p[1] = sin(params[1]) * sin(params[0]);
			p[2] = cos(params[1]);
				
			double angDist2 = ComputeSquaredAngleSum(p, 1);	
			MeasureType cost = angDist2;
			
			return cost;
		};
		
		MeasureType GetValueWithCartesianCoordinate(const ParametersType& p) const {
			double mm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
			
			double angDist2 = ComputeSquaredAngleSum(p, mm);
			MeasureType cost = angDist2 + p[3] * (1 - mm);
			
			return cost;
		};
};
};
#endif
