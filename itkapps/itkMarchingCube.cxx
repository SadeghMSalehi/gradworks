#include "itkImageCommon.h"
#include "itkBinaryMask3DSource.h"

itk::Image<unsigned int, 3> ImageType;

using namespace std;

int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << argv[0] << ": input-image output-meta" << endl;
		return;
	}
}
