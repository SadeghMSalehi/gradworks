#include "Eigen/Sparse"
#include "unsupported/Eigen/SparseExtra"
#include "iostream"

class SparseMatrix {
	protected:
		int _row;
		int _col;

	public:
		SparseMatrix(int r, int c) {
			_row = r;
			_col = c;
		}

		int Row() { return _row; }
		int Col() { return _col; }
		virtual void resize(int r, int c) = 0;
		virtual void set(int r, int c, double v) = 0;
		virtual double get(int r, int c) = 0;
		virtual void clear() = 0;
};


class EigenSparseMatrix : public SparseMatrix {
	public: 
		typedef Eigen::SparseMatrix<double, Eigen::RowMajor> EigenMatrixType;

	private:
		EigenMatrixType A;

	public:

		EigenSparseMatrix() : SparseMatrix(1,1) {};

		EigenSparseMatrix(int r, int c) : SparseMatrix(r,c) {
			A.resize(r, c);
		};

		void resize(int r, int c) {
			A.resize(r, c);
		}

		void set(int r, int c, double v) {
			A.insert(r, c) = v;
		}

		double get(int r, int c) {
			return A.coeff(r, c);
		}

		void reserve(int n) {
			A.reserve(n);
		}

		void clear() {
			A.setZero();
		}

		int count() {
			return A.nonZeros();
		}

		bool solve(Eigen::VectorXd b, Eigen::VectorXd &x) {
			Eigen::SparseLLT<EigenMatrixType> llt(A);

			if (llt.succeeded()) {
				x = llt.solve(b);
				return true;
			} else {
				std::cout << "LLT decomposition failed" << std::endl;
			}

			return false;
		}

		void print() {
			std::cout << A << std::endl;
		}

		/*
		EigenMatrixType multiply(Eigen::VectorXd &r) {
			EigenMatrixType result = A * EigenMatrixType(r);
		}
		*/
};
