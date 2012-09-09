#include "itkImageCommon.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "fstream"
#include "iostream"
#include "queue"

using namespace std;

typedef itk::Image<double,3> ImageType;
typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;

int main(int argc, char* argv[]) {
	ifstream f(argv[1]);	
	double d;
	queue<double> q;
	int x,y,z;

	if (argc < 6) {
		cout << "usage: " << argv[0] << " text_file out_image_file size_x size_y size_c" << endl;
		return 0;
	}

	x = atoi(argv[3]);
	y = atoi(argv[4]);
	z = atoi(argv[5]);

	while (true) {
		f >> d;
		if (f.eof()) {
			break;
		}
		q.push(d);
	}
	f.close();

	if (q.size() < (unsigned int) (x * y * z)) {
		cout << "Error: the number of values are less than the number of image voxels (" << x * y * z << "). " << endl;
	}

	cout << "The number of input values : " << q.size() << endl;
	cout.flush();

	ImageType::Pointer image = NewImageT<ImageType>(x,y,z);
	IteratorType iter(image, image->GetRequestedRegion());

	cout << "Creating image with size: " << image->GetLargestPossibleRegion().GetSize() << endl; 

	for (;!iter.IsAtEnd(); ++iter) {
		iter.Set(q.front());
		q.pop();
	}

	WriteImageT<ImageType>(argv[2], image);
}
