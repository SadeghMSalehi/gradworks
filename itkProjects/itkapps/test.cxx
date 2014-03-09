#include "iostream"
#include "Eigen/Dense"

using namespace std;

int main(int argc, char* argv[]) {
	Eigen::MatrixXf mat(2,2);
	mat << 1,2,3,4;
	Eigen::VectorXf vec(2);
	vec << 1,1;
	cout << mat.fullPivLu().solve(vec) << endl;
}
