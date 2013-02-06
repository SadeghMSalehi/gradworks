#include "itkImageCommon.h"
#include "itkRigid3DTransform.h"
#include "itkAffineTransform.h"
#include "itkTransformFileWriter.h"

typedef itk::AffineTransform<double> TransformType;

int main(int argc, char* argv[]) {
	TransformType::Pointer txf = TransformType::New();
	TransformType::MatrixType mat;
	
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			mat(i,j) = i*4+j;
		}
	}

	TransformType::InputPointType p;
	p[0] = 100;
	p[1] = 200;
	p[2] = 300;

	TransformType::OutputVectorType t;
	t[0] = -100;
	t[1] = -200;
	t[2] = -300;

	txf->SetCenter(p);
	txf->SetMatrix(mat);
	txf->SetTranslation(t);

	itk::TransformFileWriter::Pointer tfmW = itk::TransformFileWriter::New();
	tfmW->AddTransform(txf);
	tfmW->SetFileName(argv[1]);
	tfmW->Update();
}
