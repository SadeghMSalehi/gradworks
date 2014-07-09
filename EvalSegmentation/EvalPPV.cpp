#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkPoint.h>
// for lesion match purpose
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>

using namespace std;

typedef itk::Image<unsigned short, 3>                      SegmentationType;
typedef SegmentationType::Pointer                         SegmentationPointer;
typedef itk::ImageFileReader< SegmentationType >          SegmentationReaderType;
typedef itk::ImageFileWriter< SegmentationType >          SegmentationWriterType;
typedef itk::ImageRegionIterator<SegmentationType>        IteratorType;
typedef itk::ConnectedComponentImageFilter<SegmentationType,SegmentationType> ConnectiveFilterType;
typedef itk::RelabelComponentImageFilter<SegmentationType,SegmentationType> RelabelFilterType;


typedef itk::Point<float,3> PointType;


int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    std::cout << "Usage: " << argv[0] << " segmentation reference [-o OutputName] [-t threshold]" << std::endl;
    std::cout << "Output: resultFilename, referenceFilename, sensitivity, specificity, PPV" << std::endl;
    return -1;
  }

  // initialisation
  SegmentationReaderType::Pointer resultReader = SegmentationReaderType::New();
  SegmentationReaderType::Pointer validationReader = SegmentationReaderType::New();
  SegmentationPointer resultImage = resultReader->GetOutput();
  SegmentationPointer validationImage = validationReader->GetOutput();

  char resultFilename[512];
  sprintf( resultFilename, "%s", argv[1] );
  char referenceFilename[512];
  sprintf( referenceFilename, "%s", argv[2] );

  double threshold = 0.5;
  for (int i = 3; i < argc; i++) {
    if (strcmp(argv[i], "-t") == 0 && i < (argc - 1)) {
       threshold = atof(argv[i+1]);
    }
  }

  std::cerr << "setting threhold as " << threshold << std::endl;

  // read result image:
  try {
    //cout << "loading result image ... " << endl;
    resultReader->SetFileName( resultFilename );
    resultReader->Update();
  }
  catch(...) {
    printf( "Could not load image %s\n", resultFilename );
    return -2;
  }
  // read validation image:
  try {
    //cout << "loading reference image ... " << endl;
    validationReader->SetFileName( referenceFilename );
    validationReader->Update();
  }
  catch(...) {
    printf( "Could not load image %s\n", referenceFilename );
    return -2;
  }
  // check if images have the same size
  SegmentationType::RegionType resRegion = resultImage->GetLargestPossibleRegion();
  SegmentationType::RegionType valRegion = validationImage->GetLargestPossibleRegion();
  if (resRegion.GetSize() != valRegion.GetSize())
  {
    printf( "Image sizes for %s and %s are different!\n", resultFilename, referenceFilename );
    return -3;
  }

  SegmentationType::SpacingType valSpacing = validationImage->GetSpacing();
  double volumeFactor = 0.001*valSpacing[0]*valSpacing[1]*valSpacing[2];
  SegmentationType::SpacingType resSpacing = resultImage->GetSpacing();
  if (resSpacing[0]!=valSpacing[0] || resSpacing[1]!=valSpacing[1] || resSpacing[2]!=valSpacing[2]) 
  {
    std::cerr << "WARNING: Spacing of segmentation different from reference! Sum "
            << fabs(resSpacing[0] - valSpacing[0]) + fabs(resSpacing[1] - valSpacing[1]) + fabs(resSpacing[2] - valSpacing[2])
            << std::endl;

  }


  // sens/spec/ppv calculation
  unsigned long tp = 0, tn = 0, fp = 0, fn = 0;
  /*
  IteratorType resIt( resultImage, resRegion );
  IteratorType valIt( validationImage, valRegion );

  for ( resIt.GoToBegin(), valIt.GoToBegin(); !resIt.IsAtEnd(); ++resIt, ++valIt ) {
    if (resIt.Get()!=0) {
      if (valIt.Get()>threshold) {
		tp ++;
      } else {
		fp ++;
	  }	
    } else {
      if (valIt.Get()>threshold) {
		fn ++;
	  } else {
		tn ++;
      }
    }
  }
  */

  // Lesion comparison
  //std::cout << "Connectivity processing ..." << endl;
  ConnectiveFilterType::Pointer connectiveInputFilter = ConnectiveFilterType::New();
  ConnectiveFilterType::Pointer connectiveValidFilter = ConnectiveFilterType::New();
  RelabelFilterType::Pointer relabelInputFilter = RelabelFilterType::New();
  RelabelFilterType::Pointer relabelValidFilter = RelabelFilterType::New();

  connectiveInputFilter->SetInput(resultImage);  
  connectiveValidFilter->SetInput(validationImage);

  relabelInputFilter->SetInput(connectiveInputFilter->GetOutput());
  relabelValidFilter->SetInput(connectiveValidFilter->GetOutput());  

  relabelInputFilter->Update();
  
  SegmentationType::Pointer relabelInputImage = relabelInputFilter->GetOutput();
  SegmentationType::Pointer relabelValidImage = relabelValidFilter->GetOutput();

  /*
  SegmentationWriterType::Pointer writer = SegmentationWriterType::New();
  writer->SetInput(relabelInputImage);
  writer->SetFileName("relabelInputImage.gipl.gz");
  writer->Write();

  SegmentationWriterType::Pointer writer2 = SegmentationWriterType::New();
  writer2->SetInput(relabelValidImage);
  writer2->SetFileName("relabelValidImage.gipl.gz");
  writer2->Write();
  //std::cout << "Connectivity processing finished ..." << endl;
  */


  // lesion comparison
  IteratorType resLabelIt( relabelInputImage, resRegion ), valLabelIt( relabelValidImage, valRegion );

  int maxResultLabel = 0, maxValidLabel = 0;
  for ( resLabelIt.GoToBegin(), valLabelIt.GoToBegin(); !resLabelIt.IsAtEnd(); ++resLabelIt, ++valLabelIt ) {
	if (resLabelIt.Get() > maxResultLabel) {
	  maxResultLabel = resLabelIt.Get(); 
    }

	if (valLabelIt.Get() > maxValidLabel) {
	  maxValidLabel = valLabelIt.Get(); 
    }
  }

  maxResultLabel++;
  maxValidLabel++;

  cout << maxResultLabel << " " << maxValidLabel << endl;

  int* resultLabels = new int[maxResultLabel];
  int* validLabels  = new int[maxValidLabel];

  memset(resultLabels, 0, sizeof(int) * (maxResultLabel));
  memset(validLabels, 0, sizeof(int) * (maxValidLabel));

 for ( resLabelIt.GoToBegin(), valLabelIt.GoToBegin(); !resLabelIt.IsAtEnd(); ++resLabelIt, ++valLabelIt ) {
	// if result is lesion
    if (resLabelIt.Get() != 0) {
	  // if valid is lesion
      if (valLabelIt.Get() != 0) {
		// then true positive
        resultLabels[resLabelIt.Get()] = 2;
		validLabels[valLabelIt.Get()]  = 2;
      } else {
		// else false positive
        if (resultLabels[resLabelIt.Get()] == 0) resultLabels[resLabelIt.Get()] = 1;
	  }
    } else {
	  // if result is not lesion
      if (valLabelIt.Get() != 0) {
		// else false negative
        if (validLabels[valLabelIt.Get()] == 0) validLabels[valLabelIt.Get()] = 1;
	  } else {
		// else true negative
      }
    }
  }

  int numTP = 0, numFP = 0, numFN = 0;
  for (int i = 0; i < maxResultLabel; i++) {
 	switch (resultLabels[i]) {
		case 1:
			numFP++;
			break;
    }
  }

  for (int i = 0; i < maxValidLabel; i++) {
 	switch (validLabels[i]) {
		case 1:
			numFN++;
			break;
		case 2:
			numTP++;
			break;
    }
  }

  double truePositives  = (numTP / double(maxValidLabel - 1)) * 100;
  double falsePositives = (numFP / double(maxResultLabel - 1)) * 100;
  double falseNegatives = (numFN / double(maxValidLabel - 1)) * 100;

  if (maxResultLabel == 1) { // no segmentation in result file
	falsePositives = 0;
  }
  if (maxValidLabel == 1) { // no segmentation in reference file
	truePositives = 0;
	falseNegatives = 0;
  }

  delete[] resultLabels, validLabels; 
  
  double sensitivity = double(tp) / double(tp + fn);
  double specificity = double(tn) / double(tn + fp);
  double ppv = double(tp) / double(tp + fp);

  char resultBuffer[1024];
  sprintf(resultBuffer, "%s; %s; %f; %f; %f; %d; %d; %d; %d; %.4f; %.4f; %.4f;\n", resultFilename, referenceFilename, truePositives, falsePositives, falseNegatives, tp, fp, tn, fn, sensitivity, specificity, ppv);
  printf(resultBuffer);

  return 0;
}
