#include "itkImageIO.h"
#include "itkContinuousIndex.h"
#include "itkTransformFileWriter.h"
#include "itkTranslationTransform.h"
#include "popt.h"

typedef itk::Image<int,3> ImageType;

static const char* transformOutput;

static struct poptOption optionsTable[] = {
    { "transformOutput", (char) 't', POPT_ARG_STRING, (void*) &transformOutput, 0,
        "filename for output transformation", NULL  },
    POPT_AUTOALIAS
    POPT_AUTOHELP
    POPT_TABLEEND
};

int main(int argc, const char* argv[]) {
  poptContext ctx = poptGetContext("itkImageMatch", argc, argv, optionsTable, 0);
  int option = poptGetNextOpt(ctx);
  const char** args = poptGetArgs(ctx);
  argc = 0;
  while (args[argc++] != NULL);
  if (argc < 2) {
      poptPrintUsage(ctx, stdout, 0);
      return 0;
  }

  const char* inputImage = args[0];
  const char* referenceImage = args[1];

  itkImageIO<ImageType> io;
  ImageType::Pointer in = io.ReadImageT(inputImage);
  ImageType::Pointer out = io.ReadImageT(referenceImage);

  ImageType::RegionType inRegion = in->GetLargestPossibleRegion();
  ImageType::RegionType outRegion = out->GetLargestPossibleRegion();

  ImageType::SizeType inSz = inRegion.GetSize();
  ImageType::SizeType outSz = outRegion.GetSize();

  typedef itk::TranslationTransform<double, ImageType::ImageDimension> TransformType;
  TransformType::OutputVectorType translation;
  for (int i = 0; i < ImageType::ImageDimension; i++) {
    translation[i] = outSz[i] / 2.0 - inSz[i] / 2.0;
  }
  TransformType::Pointer transform = TransformType::New();
  transform->Translate(translation);
  transform->Update();
}
