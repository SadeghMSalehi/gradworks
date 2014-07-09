#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkPoint.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryThresholdImageFilter.h> 
#include <itkSubtractImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <ANN/ANN.h>
#include "itkMinimumMaximumImageCalculator.h"

using namespace std;



typedef itk::Image<unsigned short, 3>                      SegmentationType;
typedef SegmentationType::Pointer                         SegmentationPointer;
typedef itk::ImageFileReader< SegmentationType >          SegmentationReaderType;
typedef itk::ImageFileWriter< SegmentationType >          SegmentationWriterType;
typedef itk::ImageRegionIterator<SegmentationType>        IteratorType;
typedef itk::BinaryBallStructuringElement<char, 3>        StructuringType;
typedef itk::BinaryErodeImageFilter<SegmentationType, SegmentationType, StructuringType> ErodeFilterType;
typedef itk::BinaryThresholdImageFilter< SegmentationType, SegmentationType> ThreshFilterType;
typedef itk::SubtractImageFilter<SegmentationType, SegmentationType, SegmentationType> SubFilterType;
typedef itk::Point<float,3> PointType;

typedef itk::MinimumMaximumImageCalculator< SegmentationType > MinMaxCalcType;

// for lesion match purpose
typedef itk::ConnectedComponentImageFilter<SegmentationType,SegmentationType> ConnectiveFilterType;
typedef itk::RelabelComponentImageFilter<SegmentationType,SegmentationType> RelabelFilterType;

static int refPerc = 25;
// reference error for comparison => if a score of 90 = rater performance => refPerc 10

int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    std::cout << "Usage: " << argv[0] << " segmentation reference [-c voe rvd asd tp fp] [-o OutputName] [-r refPerfScore}" << std::endl;
    std::cout << "     -c    reference numbers, voe=Volume overlap error, rvd = relative absolute diffference percentage, asd = average surface distance,  tp = true positives, fp = false positive" << std::endl;
    std::cout << "Output: resultFilename, referenceFilename, absVolumeDifPerc, avgDistance, tanimotoError, true positives, false positives" << std::endl;
    std::cout << "     -r refPerfScore reference error for comparison => if a score of 90 is equal to a rater performance then this should be 100-90=10" << std::endl;
    return -1;
  }

  // initialisation
  SegmentationReaderType::Pointer resultReader = SegmentationReaderType::New();
  SegmentationReaderType::Pointer validationReader = SegmentationReaderType::New();
  SegmentationPointer resultImage = resultReader->GetOutput();
  SegmentationPointer validationImage = validationReader->GetOutput();
  StructuringType structuringBall;
  structuringBall.SetRadius( 1 );
  structuringBall.CreateStructuringElement();
  ThreshFilterType::Pointer threshFilter1 = ThreshFilterType::New();
  threshFilter1->SetLowerThreshold(1);
  threshFilter1->SetUpperThreshold(255);
  threshFilter1->SetOutsideValue (0);
  threshFilter1->SetInsideValue (1);
  ErodeFilterType::Pointer erodeFilter1 = ErodeFilterType::New();
  erodeFilter1->SetInput( threshFilter1->GetOutput() );
  erodeFilter1->SetKernel( structuringBall );
  erodeFilter1->SetErodeValue( 1 );
  SubFilterType::Pointer subFilter1 = SubFilterType::New();
  subFilter1->SetInput2( erodeFilter1->GetOutput() );
  ThreshFilterType::Pointer threshFilter2 = ThreshFilterType::New();
  threshFilter2->SetLowerThreshold(1);
  threshFilter2->SetUpperThreshold(255);
  threshFilter2->SetOutsideValue (0);
  threshFilter2->SetInsideValue (1);
  ErodeFilterType::Pointer erodeFilter2 = ErodeFilterType::New();
  erodeFilter2->SetInput( threshFilter2->GetOutput() );
  erodeFilter2->SetKernel( structuringBall );
  erodeFilter2->SetErodeValue( 1 );
  SubFilterType::Pointer subFilter2 = SubFilterType::New();
  subFilter2->SetInput2( erodeFilter2->GetOutput() );
  

  char resultFilename[512];
  sprintf( resultFilename, "%s", argv[1] );
  char referenceFilename[512];
  sprintf( referenceFilename, "%s", argv[2] );

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
  PointType pnt;

  // compute border pixels and init kd-tree for image 1
  threshFilter1->SetInput( resultImage );
  subFilter1->SetInput1( threshFilter1->GetOutput() );
  subFilter1->UpdateLargestPossibleRegion();
  SegmentationPointer borderImg1 = subFilter1->GetOutput();

  // compute border pixels and init kd-tree for image 2
  threshFilter2->SetInput( validationImage );
  subFilter2->SetInput1( threshFilter2->GetOutput() );
  subFilter2->UpdateLargestPossibleRegion();
  SegmentationPointer borderImg2 = subFilter2->GetOutput();

  // check if file is empty, if so create bogus single voxel at (1,1,1)
  MinMaxCalcType::Pointer minmaxCalc = MinMaxCalcType::New();
  minmaxCalc->SetImage(borderImg1);
  minmaxCalc->Compute();
  double max1 = minmaxCalc->GetMaximum();
  minmaxCalc->SetImage(borderImg2);
  minmaxCalc->Compute();
  double max2 = minmaxCalc->GetMaximum();
  double avgDistance = 0;
  double avgSqrDistance = 0;
  double maxDistance = 0;
  //std::cerr << "Max:"  << max << std::endl;
  if (max1 == 0 || max2 == 0) {
    if (max1 != 0 || max2 != 0) {
      std::cerr << "WARNING: Empty data in only result or reference image, bogus numbers for surface errors" << std::endl;
      avgDistance = valSpacing[0] * valRegion.GetSize()[0] / 2;
      avgSqrDistance = avgDistance;
      maxDistance = avgDistance;
    }
  } else {

    IteratorType it1( borderImg1, borderImg1->GetLargestPossibleRegion() );
    unsigned int numBorderPts1 = 0;
    for ( it1.GoToBegin(); !it1.IsAtEnd(); ++it1 ) {
      if (it1.Get() != 0) numBorderPts1++;
    }
    ANNpointArray borderPts1 = annAllocPts( numBorderPts1, 3 );
    numBorderPts1 = 0;
    for ( it1.GoToBegin(); !it1.IsAtEnd(); ++it1 ) {
      if (it1.Get() != 0) {
	validationImage->TransformIndexToPhysicalPoint( it1.GetIndex(), pnt );
	for (int d=0; d<3; d++) borderPts1[numBorderPts1][d] = pnt[d];
	numBorderPts1++;
      }
    }
    ANNkd_tree *borderTree1 = new ANNkd_tree( borderPts1, numBorderPts1, 3 );

    //cout << "Calculating surface distance measure.." << endl;
    IteratorType it2( borderImg2, borderImg2->GetLargestPossibleRegion() );
    unsigned int numBorderPts2 = 0;
    for ( it2.GoToBegin(); !it2.IsAtEnd(); ++it2 ) {
      if (it2.Get() != 0) numBorderPts2++;
    }
    ANNpointArray borderPts2 = annAllocPts( numBorderPts2, 3 );
    numBorderPts2 = 0;
    for ( it2.GoToBegin(); !it2.IsAtEnd(); ++it2 ) {
      if (it2.Get() != 0) {
	validationImage->TransformIndexToPhysicalPoint( it2.GetIndex(), pnt );
	for (int d=0; d<3; d++) borderPts2[numBorderPts2][d] = pnt[d];
	numBorderPts2++;
      }
    }
    ANNkd_tree *borderTree2 = new ANNkd_tree( borderPts2, numBorderPts2, 3 );
    
    // calculate surface distance measures
    ANNidxArray  nnIdx = new ANNidx[1];
    ANNdistArray dists = new ANNdist[1];
  
    for(unsigned int idx1=0; idx1<numBorderPts1; idx1++) {
      borderTree2->annkSearch( borderPts1[idx1], 1, nnIdx, dists);
      avgSqrDistance += dists[0];
      double d = sqrt( dists[0] );
      avgDistance += d;
      if (d>maxDistance) maxDistance = d;
    }
    
    for(unsigned int idx2=0; idx2<numBorderPts2; idx2++) {
      borderTree1->annkSearch( borderPts2[idx2], 1, nnIdx, dists);
      avgSqrDistance += dists[0];
      double d = sqrt( dists[0] );
      avgDistance += d;
      if (d>maxDistance) maxDistance = d;
    }

    double numBorderPts = numBorderPts1 + numBorderPts2;
    avgDistance /= numBorderPts;
    avgSqrDistance /= numBorderPts;
    avgSqrDistance = sqrt( avgSqrDistance );

    // clean up
    annDeallocPts( borderPts1 );
    annDeallocPts( borderPts2 );
    delete borderTree1;
    delete borderTree2;
    delete[] nnIdx;
    delete[] dists;
  }
  //cout << "Calculating Tanimoto overlap .." << endl;

  char resultBuffer[1024];
  // Tanimoto overlap metric
  unsigned long volume1=0, volume2=0, intersection=0;
  IteratorType resIt( resultImage, resRegion ), valIt( validationImage, valRegion );
  for ( resIt.GoToBegin(), valIt.GoToBegin(); !resIt.IsAtEnd(); ++resIt, ++valIt ) {
    if (resIt.Get()!=0) {
      volume1++;
      if (valIt.Get()!=0) {
        volume2++;
        intersection++;
      }
    }
    else {
      if (valIt.Get()!=0) volume2++;
    }
  }


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
  SegmentationWriterType::Pointer writer = SegmentationWriterType::New();
  writer->SetInput(relabelInputImage);
  writer->SetFileName("relabelInputImage.gipl.gz");
  writer->Write();

  SegmentationType::Pointer relabelValidImage = relabelValidFilter->GetOutput();
  SegmentationWriterType::Pointer writer2 = SegmentationWriterType::New();
  writer2->SetInput(relabelValidImage);
  writer2->SetFileName("relabelValidImage.gipl.gz");
  writer2->Write();
  //std::cout << "Connectivity processing finished ..." << endl;


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
  

  double tanimotoVal = 100.0 * (double)(intersection) / ((double)(volume1+volume2-intersection));
  double tanimotoError = 100.0 - tanimotoVal;
  double volumeSeg = (double)volume1 * volumeFactor;
  double volumeRef = (double)volume2 * volumeFactor;
  double volumeDif = volumeSeg - volumeRef;
  double volumeDifPerc = 100.0 * volumeDif / volumeRef;

  sprintf( resultBuffer, "%s; %s; %.4f; %.4f; %.4f; %.4f; %.4f;\n", 
	   resultFilename, referenceFilename,
	   fabs(volumeDifPerc), avgDistance, tanimotoError, truePositives, falsePositives );

  printf( resultBuffer );


  /*
  // check if metrics should be scored according to given references
  bool calcScores = false;
  double refVOE=0, refRVD=0, refASD=0, refRMSSD=0, refMSD=0, refTP=0, refFP=0;
  for (int i=3; i<argc; i++) {
    if ( strcmp( argv[i], "-c" )==0 && i<(argc-1)) {
      refVOE = atof( argv[i+1] );
      refRVD = atof( argv[i+2] );
      refASD= atof( argv[i+3] );
      refTP = atof( argv[i+4] );
	  refFP = atof( argv[i+5] );

	  // printf("Ref scores: %f %f %f %f %f\n", refVOE, refRVD, refASD, refTP, refFP);

      if (refVOE!=0 && refRVD!=0 && refASD!=0 && refTP!=0 && refFP != 0) calcScores = true;
      else std::cout << "Reference values have to be >0! Not calculating scores..." << std::endl;
      break;
    }
  }
  // reference performance
  for (int i=3; i<argc; i++) {
    if ( strcmp( argv[i], "-r" )==0 && i<(argc-1)) {
      refPerc = atoi( argv[i+1] );
      break;
    }
  }



  // save evaluation data with this name
  std::string evaluationName = "";
  for (int i=3; i<argc; i++) {
    if ( strcmp( argv[i], "-o" )==0 && i<(argc-1)) {
      evaluationName = argv[i+1];
      break;
    }
  }
  if (evaluationName == "") {
    evaluationName = "evaluation.txt";
  }
  
  std::string scoreName = "";
  if (calcScores) {
    scoreName = "scoring-" + evaluationName;
  }

   // append info to specified text file:
  FILE *file = fopen( evaluationName.c_str(), "a" );
  fputs( resultBuffer, file );
  fclose( file );

  // check if score file has to be written
  if (calcScores) {
    double scoreVOE = 100 - refPerc *( tanimotoError / refVOE );
    if (scoreVOE < 0) scoreVOE = 0;
    double scoreRVD = 100 - refPerc * ( fabs( volumeDifPerc ) / refRVD );
    if (scoreRVD < 0) scoreRVD = 0;
    double scoreASD = 100 - refPerc * ( avgDistance / refASD );
    if (scoreASD < 0) scoreASD = 0;
    double scoreTruePositives = refPerc * (truePositives / refTP);
	if (scoreTruePositives < 0) scoreTruePositives = 0;
    double scoreFalsePositives = 100 - refPerc * (falsePositives / refFP);
	if (scoreFalsePositives < 0) scoreFalsePositives = 0;

    double totalScore = (scoreVOE+scoreRVD+scoreASD+scoreTruePositives+scoreFalsePositives) / 5;
    sprintf( resultBuffer, "%s; %.f; %.f; %.f; %.f; %.f; %.f;\n", 
      resultFilename, scoreVOE, scoreRVD, scoreASD, scoreTruePositives, scoreFalsePositives, totalScore );
    file = fopen( scoreName.c_str(), "a" );
    fputs( resultBuffer, file );
    fclose( file );
  }
  */


  return 0;
}
