#include "itkImage.h"
#include "itkImageCommon.h"
#include "vector"
#include "itkImageRegionIteratorWithIndex.h"
#include "iostream"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "boost/program_options.hpp"
#include "vtkFloatArray.h"
#include "vtkIntArray.h"

using namespace std;
namespace po = boost::program_options;

typedef itk::Vector<double,3> VectorType;
typedef itk::Image<VectorType,3> ImageType;

typedef itk::Image<short,3> LabelImageType;
typedef itk::ImageRegionIteratorWithIndex<LabelImageType> LabelIteratorType;

void extractNeighorVectors(string &infile, vector<float> &xyz, string &outfile) {
	int ret = 0;	
	ImageType::Pointer p = ReadImageT<ImageType>(infile.c_str(), ret);

	ImageType::IndexType idx;
	idx[0] = int(xyz[0]);
	idx[1] = int(xyz[1]);
	idx[2] = int(xyz[2]);

	vtkPoints* pts = vtkPoints::New();
	ImageType::IndexType nbr;
	for (int i = -1; i < 2; i++) {
		nbr[2] = idx[2] + i;	
		for (int j = -1; j < 2; j++) {
			nbr[1] = idx[1] + j;	
			for (int k = -1; k < 2; k++) {
				nbr[0] = idx[0] + k;	
				VectorType v = p->GetPixel(nbr);
				pts->InsertNextPoint(v[0], v[1], v[2]);
			}
		}
	}

	vtkPolyData* pdata = vtkPolyData::New();
	pdata->SetPoints(pts);

	vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
	writer->SetFileName(outfile.c_str());
	writer->SetInput(pdata);
	writer->Write();
}

void extractSegmentationVectors(string &infile, string &segfile, string &outfile, bool isVectorField) {
	int ret = 0;	
	ImageType::Pointer p = ReadImageT<ImageType>(infile.c_str(), ret);

	LabelImageType::Pointer sImg = ReadImageT<LabelImageType>(segfile.c_str(), ret);
	LabelIteratorType iter(sImg, sImg->GetRequestedRegion());

	vtkIntArray* labels = vtkIntArray::New();
	vtkFloatArray* vectors = vtkFloatArray::New();

	labels->SetNumberOfComponents(1);
	labels->SetName("Labels");

	vectors->SetName("Eigenvectors");
	vectors->SetNumberOfComponents(3);

	vtkPoints* pts = vtkPoints::New();
	for (; !iter.IsAtEnd(); ++iter) {
		short l = iter.Get();
		if (l > 0) {
			LabelImageType::IndexType nbr = iter.GetIndex();
			VectorType v = p->GetPixel(nbr);
			if (!isVectorField) {
				pts->InsertNextPoint(v[0], v[1], v[2]);
				labels->InsertNextValue(int(l));
			} else {
				pts->InsertNextPoint(nbr[0], nbr[1], nbr[2]);
				labels->InsertNextValue(int(l));
				vectors->InsertNextTuple3(v[0], v[1], v[2]);
			}
		}
	}

	vtkPolyData* pdata = vtkPolyData::New();
	pdata->SetPoints(pts);
	pdata->GetPointData()->SetScalars(labels);

	if (isVectorField) {
		pdata->GetPointData()->SetVectors(vectors);
	}

	vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
	writer->SetFileName(outfile.c_str());
	writer->SetInput(pdata);
	writer->Write();
}
	
int main(int argc, char* argv[]) {
	string infile, outfile, segfile;
	vector<float> xyzCoord;

	po::options_description desc("Options");
	desc.add_options()
		("help", "print help messages")
		("xyz", po::value<vector<float> >(&xyzCoord), "x,y,z coordinate to look up neighbor eigen vectors")
		("segmentation", po::value<std::string>(&segfile), "segmentation file to extract eigen vectors")
		("vectorfield", "generate vector field instead of orientation map")
	;

	po::options_description args("Args");
	args.add_options()
		("input", po::value<string>(&infile), "input image")
		("output", po::value<string>(&outfile), "output image")
	;

	po::options_description allOptions("All Options");
	allOptions.add(args).add(desc);

	po::positional_options_description p;
	p.add("input", 1);
	p.add("output", 1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(allOptions).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("input") < 1 || vm.count("output") < 1) {
		cout << "usage: " << argv[0] << " [options] output-vti" << endl;
		cout << desc << endl;
		return 0;
	}

	if (vm.count("xyz") > 0) {
		extractNeighorVectors(infile, xyzCoord, outfile);
	} else if (vm.count("segmentation") > 0) {
		extractSegmentationVectors(infile, segfile, outfile, vm.count("vectorfield") > 0);
	}
}
