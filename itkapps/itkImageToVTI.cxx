#include "itkImageCommon.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkStructuredGrid.h"
#include "vtkStructuredGridWriter.h"
#include "vtkImageData.h"
#include "vtkDoubleArray.h"
#include "vtkXMLImageDataWriter.h"
#include "boost/program_options.hpp"

#include "iostream"


using namespace std;
namespace po = boost::program_options;

typedef itk::Image<unsigned short,3> LabelImageType;
typedef itk::Vector<double,3> VectorType;
typedef itk::Image<VectorType,3> VectorImageType;
typedef itk::Image<double,3> DoubleImageType;
typedef itk::ImageRegionIteratorWithIndex<VectorImageType> IteratorType;
typedef itk::ImageRegionIteratorWithIndex<LabelImageType> LabelIteratorType;

void SetArrayTuple(vtkDataArray* a, int i, VectorType v) {
	for (int j = 0; j < a->GetNumberOfComponents(); j++) {
		a->SetComponent(i, j, v[j]);
	}
}

void SetArrayTuple(vtkDataArray* a, int i, double v) {
	a->SetTuple1(i, v);
}

template <class X>
void ConvertImageT(string& imageFile, vtkImageData* imgData, const char* attrName, int numberOfComponents) {
	int ret = 0;
  typename X::Pointer srcImg = ReadImageT<X>(imageFile.c_str(), ret);
	typename X::SizeType srcSize = srcImg->GetRequestedRegion().GetSize();
	typename X::PointType srcOrigin = srcImg->GetOrigin();
	typename X::SpacingType srcSpacing = srcImg->GetSpacing();

	imgData->SetOrigin(srcOrigin[0], srcOrigin[1], srcOrigin[2]);
	imgData->SetSpacing(srcSpacing[0], srcSpacing[1], srcSpacing[2]);
	imgData->SetDimensions(srcSize[0], srcSize[1], srcSize[2]);

	const int nPoints = srcSize[0] * srcSize[1] * srcSize[2];

	vtkDoubleArray* attr = vtkDoubleArray::New();
	attr->SetNumberOfComponents(numberOfComponents);
	attr->SetName(attrName);
	attr->SetNumberOfTuples(nPoints);

	vtkPointData* pdata = imgData->GetPointData();

	switch (numberOfComponents) {
		case 1:
			pdata->SetScalars(attr);
			break;
		case 3:
			pdata->SetVectors(attr);
			break;
		case 9:
			pdata->SetTensors(attr);
			break;
		default:
			pdata->AddArray(attr);
			break;
	}

	int cnt = 0;
	for (unsigned int z = 0; z < srcSize[2]; z++) {
		for (unsigned int y = 0; y < srcSize[1]; y++) {
			for (unsigned int x = 0; x < srcSize[0]; x++) {
				typename X::IndexType idx;
				idx[0] = x;
				idx[1] = y;
				idx[2] = z;
				typename X::PixelType v = srcImg->GetPixel(idx);
				SetArrayTuple(attr, cnt, v);
				cnt ++;
			}
		}
	}
}

int main(int argc, char* argv[]) {
	string vectorImage, scalarImage, outputImage;

	po::options_description desc("Options");
	desc.add_options()
		("help", "print help messages")
		("vector,v", po::value<string>(&vectorImage), "Vector attribute image")
		("scalar,s", po::value<string>(&scalarImage), "Scalar attribute image")
	;

	po::options_description args("Args");
	args.add_options()
		("output-vti", po::value<string>(&outputImage), "output VTI image")
	;

	po::options_description allOptions("All Options");
	allOptions.add(args).add(desc);

	po::positional_options_description p;
	p.add("output-vti", 1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(allOptions).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("output-vti") < 1 || (vm.count("vector") < 1 && vm.count("scalar") < 1)) {
		cout << "usage: " << argv[0] << " [options] output-vti" << endl;
		cout << desc << endl;
		return 0;
	}

	vtkImageData* imgData = vtkImageData::New();

	if (vectorImage != "") {
		ConvertImageT<VectorImageType>(vectorImage, imgData, "vector", 3);
		cout << "Set vector: " << scalarImage << endl;
	}
	if (scalarImage != "") {
		ConvertImageT<DoubleImageType>(scalarImage, imgData, "scalar", 1);
		cout << "Set scalar: " << scalarImage << endl;
	}

	vtkXMLImageDataWriter* w = vtkXMLImageDataWriter::New();
	w->SetFileName(outputImage.c_str());
	w->SetDataModeToAppended();
	w->EncodeAppendedDataOff();
	w->SetCompressorTypeToZLib();
	w->SetInput(imgData);
	w->Write();
}

/*
int oldmain(int argc, char* argv[]) {
  if (argc < 4) {
    cout << "usage: " << argv[0] << " input-image mask-image output-vti" << endl;
    return 0;
  }

  int ret = 0;
  VectorImageType::Pointer srcImg = ReadImageT<VectorImageType>(argv[1], ret);
  LabelImageType::Pointer labelImg = ReadImageT<LabelImageType>(argv[2], ret);
  
  vtkUnstructuredGrid* vtkGrid = vtkUnstructuredGrid::New();
  vtkPoints* points = vtkPoints::New();
  vtkDoubleArray* eigenVectors = vtkDoubleArray::New();
  eigenVectors->SetNumberOfComponents(3);
  eigenVectors->SetName("EigenVectors");
  
  LabelIteratorType iter(labelImg, labelImg->GetRequestedRegion());
  for (; !iter.IsAtEnd(); ++iter) {
    if (iter.Get() > 0) {
      LabelImageType::IndexType iterIdx = iter.GetIndex();
      VectorType v = srcImg->GetPixel(iterIdx);
      points->InsertNextPoint(iterIdx[0], iterIdx[1], iterIdx[2]);
      eigenVectors->InsertNextTuple3(v[0], v[1], v[2]);
    }
  }
  vtkGrid->SetPoints(points);
  vtkGrid->GetPointData()->SetVectors(eigenVectors);

  vtkUnstructuredGridWriter* gridWriter = vtkUnstructuredGridWriter::New();
  gridWriter->SetInput(vtkGrid);
  gridWriter->SetFileName(argv[3]);
  gridWriter->Write();

  return 0;
}
*/
