#include "itkTransformFileWriter.h"
#include "itkMatrixOffsetBaseTransform.h"
#include "iostream"

using namespace std;

int main(int argc, char* argv[]) {
	const char* transformIn = argv[1];
	const char* matrixOut = argv[2];

	typedef itk::TransformFileReader TransformReaderType;
	TransformReaderType::Pointer reader = TransformReaderType::New();
	reader->SetInput(transformIn);
	itk::TransformBase::Pointer transformBase = reader->GetOutput();
	typedef itk::MatrixOffsetBaseTransform<double,3,3> MatrixTramsformType;
	MatrixTramsformType* baseMatrix = dynamic_cast<MatrixTransformType*>(transformBase.GetPointer());
	MatrixTransformType::MatrixType baseMatrixRep;
	baseMatrixRep = baseMatrix->GetMatrix();
	cout << baseMatrixRep << endl;
}
