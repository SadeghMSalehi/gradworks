#include "itkImageIO.h"
#include "itkAutoCropLabelMapFilter.h"
#include "itkTransformFileWriter.h"
#include "itkSimilarity3DTransform.h"

using namespace std;
using namespace itk;
using namespace itkcmds;

typedef Image<float,3> ImageType;
typedef Similarity3DTransform<double> TransformType;

int main(int argc, char* argv[]) {
  TransformFileWriter::Pointer tfWriter = TransformFileWriter::New(); 
  TransformType::Pointer transform = TransformType::New();
  tfWriter->SetFileName("tert.txt");
  tfWriter->SetInput(transform);
  tfWriter->Update();
}
