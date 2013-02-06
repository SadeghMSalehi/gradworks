// to use multi-thread feature in itk
// create new filter class

#include "itkUnaryFunctorImageFilter.h"
#include "math.h"

#include "iostream"

using namespace std;

typedef itk::Image<double,3> ImageType;
typedef itk::Vector<double,3> VectorType; 
typedef itk::Image<VectorType,3> VectorImageType; 

namespace itk {
template<class TOutputPixel>
class VectorComponent {
	private:
		int _component;

	public:
		VectorComponent() {}
		~VectorComponent() {}

		void SetComponent(int n) {
			_component = n;
		}

		bool operator!=(const VectorComponent& other) const {
			return false;
		}

		bool operator==(const VectorComponent& other) const {
			return !(*this != other);
		}

		inline double operator()(VectorType v) const {
			return v[_component];
		}
};


class ITK_EXPORT VectorComponentFilter 
		: public UnaryFunctorImageFilter<VectorImageType,SphericalImageType,CartesianToSphericalCoordConverter> {

	public:
		typedef CartesianToSphericalCoordFilter Self;
		typedef UnaryFunctorImageFilter<VectorImageType,SphericalImageType,CartesianToSphericalCoordFilter> Superclass;
		typedef SmartPointer<Self> Pointer;
		typedef SmartPointer<const Self> ConstPointer;

		itkNewMacro(Self);
		itkTypeMacro(CartesianToSphericalCoordFilter, UnaryFunctorImageFilter);
	
	protected:
		CartesianToSphericalCoordFilter() {}
		~CartesianToSphericalCoordFilter() {}

	private:
		CartesianToSphericalCoordFilter(const Self&);
		void operator=(const Self&);

};



};
