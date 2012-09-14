#include "iostream"
#include <itkResampleImageFilter.h>
#include <itkAffineTransform.h>
#include <itkMath.h>
#include "itkMathCode.h"
#include "itkImageIO.h"
#include <itkSimilarity3DTransform.h>
#include <itkRigid3DTransform.h>
#include "itkResampler.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"

#include "itkTransformFactory.h"
#include "list"

using namespace std;
using namespace itk;
using namespace itkcmds;

typedef itk::Image<double,3> ImageType;
typedef itk::Similarity3DTransform<double> TransformType;
typedef itk::TransformFileReader TransformReaderType;
typedef itk::TransformFileWriter TransformWriterType;

int main(int argc, const char* argv[]) {
    TransformReaderType::Pointer txReader = TransformReaderType::New();
    txReader->SetFileName("txparams.txt");
    txReader->Update();
    TransformReaderType::TransformListType* txList = txReader->GetTransformList();
    TransformReaderType::TransformPointer firstTransform = txList->front();

    cout << firstTransform << endl;
    
    itkResampler<ImageType, TransformReaderType::TransformType> resampler;
    resampler.SetTransform(firstTransform);
    TransformType::ParametersType params = firstTransform->GetParameters();
    cout << firstTransform->GetNameOfClass() << " Parameters: " << params << endl;

    resampler.Run(argc, argv);
}
