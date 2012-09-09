#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include "Eigen/Sparse"
#include "unsupported/Eigen/SparseExtra"
#include "iostream"

using namespace std;
using namespace Eigen;

#define N 100

typedef SparseMatrix<double> MatrixType;

int main(int argc, char* argv[]) {

	MatrixType A(N,N);
	VectorXd b(N), x(N);

	A.reserve(N);
	for (int i = 0; i < N; i++) {
		A.insert(i,i) = (i+1)*(i+1);
		b[i] = i;
	}
	A.finalize();

	cout << b << endl;
//	cout << A << endl;

	SparseLDLT<MatrixType> ldltA(A);
	if (!ldltA.succeeded()) {
		cout << "failed" << endl;
	} 
	x = ldltA.solve(b);
	cout << x << endl;
	
}
