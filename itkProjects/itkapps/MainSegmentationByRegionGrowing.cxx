#include <vector>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <string>
#include "deque"
#include "itkImageCommon.h"
#include "DCRegionGrowingAlgorithm.h"
#include "math.h"
#include "boost/program_options.hpp"


namespace po = boost::program_options;

using namespace std;
using namespace niral;

int main(int argc, char* argv[]) {
	bool useFeature1 = false;
	bool useFeature2 = false;

  string reasonForTermination;

	vector<double> defaultAngleSearch(3,2);
	defaultAngleSearch[1] = 10;
	defaultAngleSearch[2] = 0.1;

	vector<double> defaultFeature1Values(3);
	defaultFeature1Values[0] = 0.45;
	defaultFeature1Values[1] = 1;
	defaultFeature1Values[2] = 1;

	vector<double> defaultFeature2Values(3);
	defaultFeature2Values[0] = -0.15;
	defaultFeature2Values[1] = 0.15;
	defaultFeature2Values[2] = 0.012;

	po::options_description desc("Options");
	desc.add_options()
		("help", "produce help messages")
		("roi-image", po::value<string>(), "Region of Interest [optional]")		
		("feature1-values,1", po::value<vector<double> >()
				->default_value(defaultFeature1Values, "0.45,1,1")->composing(), 
				"min-threshold,max-threshold,similarity.")
		("feature2-values,2", po::value<vector<double> >()
				->default_value(defaultFeature2Values, "-0.15,0.15,0.012")->composing(), 
				"min-threshold,max-threshold,similarity.")
		("feature1-image", po::value<string>(), "feature1 image")
		("feature2-image", po::value<string>(), "feature2 image")
		("angle-search,a", po::value<vector<double> >()
				->default_value(defaultAngleSearch, "2,10,0.1")->composing(), 
				"starting angle, upper limit, step. ie. 2,10,0.1")
    ("reason-for-termination", po::value<string>(&reasonForTermination), 
      "reason for termination; 1: angle 2: out of feature1 threshold 3: out of feature1 similarity 4: out of feature2 threshold 5: out of feature2 similarity")
		("max-pixels", po::value<int>()->default_value(100000), "maximum number of pixels for segmentation")
	;
	
	po::options_description args("Arguments");
	args.add_options()
		("input-files", po::value< vector<string> >(), "input files [input-image seed-image]")
		("output-file", po::value< string >(), "output image")
	;
	
	po::options_description all_opts("All options");
	all_opts.add(args).add(desc);
	
	po::positional_options_description p;
	p.add("input-files", 2);
  p.add("output-file", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(all_opts).positional(p).run(), vm);
	po::notify(vm);
	
	if (vm.count("input-files") < 1 || vm.count("output-file") < 1 || vm.count("help") > 0) {
		cout << "Usage: " << argv[0] << " eigen-image seed-image output-image [options]" << endl;
		cout << desc << endl;
		return 0;
	}

	int ret = 0;
	
	vector<string> inputFiles = vm["input-files"].as< vector<string> >();

	vector<double> angleSearch = vm["angle-search"].as<vector<double> >();
	vector<double> feature1Values, feature2Values;
	string feature1Image, feature2Image;

	if (vm.count("feature1-image") > 0) {
		useFeature1 = true;
		feature1Values = vm["feature1-values"].as<vector<double> >();
		feature1Image = vm["feature1-image"].as<string>();
    cout << "Feature #1: " << feature1Image << " (" << feature1Values[0] << " ~ " << feature1Values[1] << " ; " << feature1Values[2] << ")" << endl;
	}
	if (vm.count("feature2-image") > 0) {
		useFeature2 = true;
		feature2Values = vm["feature2-values"].as<vector<double> >();
		feature2Image = vm["feature2-image"].as<string>();
    cout << "Feature #2: " << feature2Image << " (" << feature2Values[0] << " ~ " << feature2Values[1] << " ; " << feature2Values[2] << ")" << endl;
	}

	string outputFile = vm["output-file"].as<string>();
	
	
	typedef itk::Image<unsigned short, 3> LabelImageType;

	VectorImageType::Pointer inputImg = ReadImageT<VectorImageType>(inputFiles[0].c_str(), ret);
	LabelImageType::Pointer seedImg  = ReadImageT<LabelImageType>(inputFiles[1].c_str(), ret);


	DCRegionGrowingAlgorithm<LabelImageType> dcrga(inputImg, seedImg);

	if (useFeature1) {
		FAImageType::Pointer fImg  = ReadImageT<FAImageType>(feature1Image.c_str(), ret);
		dcrga.SetFeature1Image(fImg);
		dcrga.SetFeature1Values(feature1Values[0], feature1Values[1], feature1Values[2]);
	}

	if (useFeature2) {
		FAImageType::Pointer fImg  = ReadImageT<FAImageType>(feature2Image.c_str(), ret);
		dcrga.SetFeature2Image(fImg);
		dcrga.SetFeature2Values(feature2Values[0], feature2Values[1], feature2Values[2]);
	}

  if (vm.count("roi-image") > 0) {
    LabelImageType::Pointer roiImg = ReadImageT<LabelImageType>(vm["roi-image"].as<string>().c_str(), ret);
    dcrga.SetROI(roiImg);
  }

	dcrga.Init();
  dcrga.SetMaxNumberOfPixels(vm["max-pixels"].as<int>());

  int outputLabel = 1;
  for (double angleInput = angleSearch[0]; angleInput < angleSearch[1]; angleInput += angleSearch[2]) {
    double threshold = cos(angleInput * M_PI / 180.0);
    cout << "Trying DCt = " << threshold << 
				" (angle = " << angleInput << " degree; outputLabel = " << outputLabel << ")" << endl; 
		dcrga.SetDC(threshold);
    dcrga.SetOutputLabel(outputLabel++);
    dcrga.DoRegionGrowing();

    cout << "# of included pixels: " << dcrga.GetVoxelCount() << endl;
    if (dcrga.IsOverPixels() || !dcrga.Reinitialize()) {
      break;
    }
  }

	if (dcrga.GetVoxelCount() > 100) {
		WriteImageT<LabelImageType>(outputFile.c_str(), dcrga.GetOutput());
	}

  if (reasonForTermination != "") {
		WriteImageT<LabelImageType>(reasonForTermination.c_str(), dcrga.GetReasonForTermination());
  }
}
