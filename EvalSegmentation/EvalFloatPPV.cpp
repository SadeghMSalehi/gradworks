#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkPoint.h>

using namespace std;



typedef itk::Image<unsigned short, 3>                      SegmentationType;
typedef SegmentationType::Pointer                         SegmentationPointer;
typedef itk::ImageFileReader< SegmentationType >          SegmentationReaderType;
typedef itk::ImageFileWriter< SegmentationType >          SegmentationWriterType;
typedef itk::ImageRegionIterator<SegmentationType>        IteratorType;

typedef itk::Image<double, 3>     		                  FloatSegmentationType;
typedef FloatSegmentationType::Pointer                         FloatSegmentationPointer;
typedef itk::ImageFileReader< FloatSegmentationType >          FloatSegmentationReaderType;
typedef itk::ImageFileWriter< FloatSegmentationType >          FloatSegmentationWriterType;
typedef itk::ImageRegionIterator<FloatSegmentationType>        FloatIteratorType;



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
  FloatSegmentationReaderType::Pointer validationReader = FloatSegmentationReaderType::New();
  SegmentationPointer resultImage = resultReader->GetOutput();
  FloatSegmentationPointer validationImage = validationReader->GetOutput();

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
  FloatSegmentationType::RegionType valRegion = validationImage->GetLargestPossibleRegion();
  if (resRegion.GetSize() != valRegion.GetSize())
  {
    printf( "Image sizes for %s and %s are different!\n", resultFilename, referenceFilename );
    return -3;
  }

  FloatSegmentationType::SpacingType valSpacing = validationImage->GetSpacing();
  double volumeFactor = 0.001*valSpacing[0]*valSpacing[1]*valSpacing[2];
  SegmentationType::SpacingType resSpacing = resultImage->GetSpacing();
  if (resSpacing[0]!=valSpacing[0] || resSpacing[1]!=valSpacing[1] || resSpacing[2]!=valSpacing[2]) 
  {
    std::cerr << "WARNING: Spacing of segmentation different from reference! Sum "
            << fabs(resSpacing[0] - valSpacing[0]) + fabs(resSpacing[1] - valSpacing[1]) + fabs(resSpacing[2] - valSpacing[2])
            << std::endl;

  }


  // Tanimoto overlap metric
  unsigned long tp = 0, tn = 0, fp = 0, fn = 0;
  IteratorType resIt( resultImage, resRegion );
  FloatIteratorType valIt( validationImage, valRegion );

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


  double sensitivity = double(tp) / double(tp + fn);
  double specificity = double(tn) / double(tn + fp);
  double ppv = double(tp) / double(tp + fp);

  char resultBuffer[1024];
  sprintf(resultBuffer, "%s; %s; %d; %d; %d; %d; %.4f; %.4f; %.4f;\n", resultFilename, referenceFilename, tp, fp, tn, fn, sensitivity, specificity, ppv);
  printf(resultBuffer);

  return 0;
}
