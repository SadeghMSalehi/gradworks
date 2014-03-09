#include "itkImageCommon.h"

typedef itk::Image<short,3> ShortImageType;
typedef itk::Image<double,3> DoubleImageType;

int main(int argc, char* argv[]) {
	ReadImageT<ShortImageType>();

}
