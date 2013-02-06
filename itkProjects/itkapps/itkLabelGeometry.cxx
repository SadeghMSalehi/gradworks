#include "itkImageCommon.h"

#include "itkLabelGeometryImageFilter.h"
#include "iostream"
#include "fstream"
#include "sstream"

typedef itk::Image<short,3> LabelMapImageType;
typedef itk::Image<unsigned short,3> ImageType;
typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
typedef itk::LabelGeometryImageFilter<LabelMapImageType,ImageType> LabelGeometryImageFilterType;

using namespace std;

int main(int argc, char* argv[]) {
	stringstream sout(stringstream::in | stringstream::out);
  int ret;

  LabelMapImageType::Pointer labelImg = ReadImageT<LabelMapImageType>(argv[1],ret);
  ImageType::Pointer inputImg = ReadImageT<ImageType>(argv[2],ret);

  LabelGeometryImageFilterType::Pointer filter = LabelGeometryImageFilterType::New();
  filter->SetInput(labelImg);
  filter->Update();

  ImageType::RegionType region = inputImg->GetRequestedRegion();
  ImageType::SizeType size = region.GetSize();
  ImageType::Pointer newImg = NewImageT<ImageType>(size[0], size[1], size[2]);
  newImg->SetSpacing(inputImg->GetSpacing());

	sout << "{" << endl;
  for (int j = 4; j < argc; j++) {
		cout << "Processing Label: " << (j-3) << endl;
    int label = atoi(argv[j]);
		if (j > 4) {
			sout << "," << endl;
		}
		sout << "\"label" << (j-3) << "\":";

		filter->CalculateOrientedBoundingBoxOn();
		filter->Update();

    LabelGeometryImageFilterType::BoundingBoxVerticesType obbList 
			= filter->GetOrientedBoundingBoxVertices(label);

    LabelGeometryImageFilterType::LabelPointType obbSize 
			= filter->GetOrientedBoundingBoxSize(label);

    LabelGeometryImageFilterType::LabelPointType obbOrigin 
			= filter->GetOrientedBoundingBoxOrigin(label);

		double labelVolume = filter->GetOrientedBoundingBoxVolume(label);

		LabelGeometryImageFilterType::MatrixType rotationMatrix
			= filter->GetRotationMatrix(label);

		LabelGeometryImageFilterType::MatrixType eigenVectors
			= filter->GetEigenvectors(label);

		LabelGeometryImageFilterType::VectorType eigenValues
			= filter->GetEigenvalues(label);

		// output as json
		sout << "{ \"obb_vertices\": {";
		for (int i = 0; i < obbList.size(); i++) {
			if (i > 0) {
				sout << ",";
			}
			sout << "\"vertice" << i << "\":[" << obbList[i][0] << "," << obbList[i][1] << "," << obbList[i][2]  << "]";
		}
		sout << "}, \"obb_origin\":[";
		for (int i = 0; i < LabelGeometryImageFilterType::LabelPointType::PointDimension; i++) {
			if (i > 0) {
				sout << ",";
			}
			sout << obbOrigin[i];
		}
		sout << "]";
		sout << ",";
		sout << "\"obb_rotation\":[";
		for (int i = 0; i < rotationMatrix.rows(); i++) {
			if (i > 0) {
				sout << ",";
			}
			for (int j = 0; j < rotationMatrix.cols(); j++) {
				if (j > 0) {
					sout << ",";
				}
				sout << rotationMatrix(i,j);
			}
		}
		sout << "]";
		sout << ",";
		sout << "\"obb_length\":[";
		for (int i = 0; i < LabelGeometryImageFilterType::LabelPointType::PointDimension; i++) {
			if (i > 0) {
				sout << ",";
			}
			sout << obbSize[i];
		}
		sout << "]";
		sout << "," ;
		sout << "\"eigenvectors\":[";
		for (int i = 0; i < eigenVectors.rows(); i++) {
			if (i > 0) {
				sout << ",";
			}
			for (int j = 0; j < eigenVectors.cols(); j++) {
				if (j > 0) {
					sout << ",";
				}
				sout << eigenVectors(i,j);
			}
		}
		sout << "]";
		sout << ",";
		sout << "\"eigenvalues\":[";
		for (int i = 0; i < eigenValues.size(); i++) {
			if (i > 0) {
				sout << ",";
			}
			sout << eigenValues[i];
		}
		sout << "]";

		sout << ",";
		sout << "\"volume\":";
		sout << labelVolume ;


		sout << "}";

//		cout << eigenVectors << endl;
//		cout << obbSize << endl;
//		cout << obbOrigin << endl;

  }
	sout << "\n}";

	ofstream jsonOutput(argv[3]);
	jsonOutput << sout.str();
	jsonOutput.close();

}
