// to use multi-thread feature in itk
// create new filter class

#include "itkUnaryFunctorImageFilter.h"
#include "math.h"

#include "iostream"

using namespace std;

typedef itk::Vector<double,2> SphericalCoordType;
typedef itk::Image<SphericalCoordType,3> SphericalImageType;
typedef itk::Vector<double,3> VectorType; 
typedef itk::Image<VectorType,3> VectorImageType; 

namespace itk {
class CartesianToSphericalCoordConverter {
	public:
		CartesianToSphericalCoordConverter() {}
		~CartesianToSphericalCoordConverter() {}

		bool operator!=(const CartesianToSphericalCoordConverter& other) const {
			return false;
		}

		bool operator==(const CartesianToSphericalCoordConverter& other) const {
			return !(*this != other);
		}

		inline SphericalCoordType operator()(VectorType v) const {
			SphericalCoordType s;
			s[0] = atan2(v[2], v[0]);//phi
			s[1] = acos(v[1]); //theta
			return s;
		}
};

class SphericalToCartesianCoordConverter {
	public:
		SphericalToCartesianCoordConverter() {}
		~SphericalToCartesianCoordConverter() {}

		bool operator!=(const SphericalToCartesianCoordConverter& other) const {
			return false;
		}

		bool operator==(const SphericalToCartesianCoordConverter& other) const {
			return !(*this != other);
		}

		inline VectorType operator()(SphericalCoordType s) const {
			VectorType v;
			double phi = s[0];
			double theta = s[1];
			v[0] = sin(theta) * cos(phi);
			v[1] = cos(theta);
			v[2] = sin(theta) * sin(phi);
			return v;
		}
};

class CartesianToHalfSphericalCoordConverter {
	public:
		CartesianToHalfSphericalCoordConverter() {}
		~CartesianToHalfSphericalCoordConverter() {}

		bool operator!=(const CartesianToHalfSphericalCoordConverter& other) const {
			return false;
		}

		bool operator==(const CartesianToHalfSphericalCoordConverter& other) const {
			return !(*this != other);
		}

		inline SphericalCoordType operator()(VectorType v) const {
			SphericalCoordType s;
			if (v[2] < 0) {
				s[0] = atan2(-v[2], -v[0]);//phi
			} else {
				s[0] = atan2(v[2], v[0]);//phi
			}
			s[1] = acos(v[1]); //theta
			return s;
		}
};


class ITK_EXPORT CartesianToSphericalCoordFilter 
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



class ITK_EXPORT CartesianToHalfSphericalCoordFilter 
		: public UnaryFunctorImageFilter<VectorImageType,SphericalImageType,CartesianToHalfSphericalCoordConverter> {

	public:
		typedef CartesianToHalfSphericalCoordFilter Self;
		typedef UnaryFunctorImageFilter<VectorImageType,SphericalImageType,CartesianToHalfSphericalCoordFilter> Superclass;
		typedef SmartPointer<Self> Pointer;
		typedef SmartPointer<const Self> ConstPointer;

		itkNewMacro(Self);
		itkTypeMacro(CartesianToHalfSphericalCoordFilter, UnaryFunctorImageFilter);
	
	protected:
		CartesianToHalfSphericalCoordFilter() {}
		~CartesianToHalfSphericalCoordFilter() {}

	private:
		CartesianToHalfSphericalCoordFilter(const Self&);
		void operator=(const Self&);

};


class ITK_EXPORT SphericalToCartesianCoordFilter 
		: public UnaryFunctorImageFilter<SphericalImageType,VectorImageType,SphericalToCartesianCoordConverter> {

	public:
		typedef SphericalToCartesianCoordFilter Self;
		typedef UnaryFunctorImageFilter<VectorImageType,SphericalImageType,SphericalToCartesianCoordFilter> Superclass;
		typedef SmartPointer<Self> Pointer;
		typedef SmartPointer<const Self> ConstPointer;

		itkNewMacro(Self);
		itkTypeMacro(SphericalToCartesianCoordFilter, UnaryFunctorImageFilter);
	
	protected:
		SphericalToCartesianCoordFilter() {}
		~SphericalToCartesianCoordFilter() {}

	private:
		SphericalToCartesianCoordFilter(const Self&);
		void operator=(const Self&);

};

};
