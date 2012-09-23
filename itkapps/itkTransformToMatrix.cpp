#include "itkTransformFileReader.h"
#include "itkMatrixOffsetTransformBase.h"
#include "iostream"

using namespace std;

int main(int argc, char* argv[]) {
	const char* transformIn = argv[1];
    
	typedef itk::TransformFileReader TransformReaderType;
	TransformReaderType::Pointer reader = TransformReaderType::New();
	reader->SetFileName(transformIn);
    reader->Update();
    TransformReaderType::TransformListType* transformList = reader->GetTransformList();
    if (transformList->size() < 1) {
        cout << "Empty transform ... " << transformList->size() << endl;
        return -2;
    }
	typedef itk::MatrixOffsetTransformBase<double,3,3> MatrixTransformType;
	MatrixTransformType* baseMatrix = dynamic_cast<MatrixTransformType*>(transformList->front().GetPointer());
    if (baseMatrix == NULL) {
        cout << "Can't convert into MatrixTransform" << endl;
        return -1;
    }
    MatrixTransformType::MatrixType baseMatrixRep = baseMatrix->GetMatrix();
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i > 0 || j > 0) {
                cout << ",";
            }
            cout << baseMatrixRep[i][j];
        }
        cout << "," << baseMatrix->GetTranslation()[i];
    }
    cout << ",0,0,0,1" << endl;
}
