/*
 * main.cpp
 *
 *  Created on: Jul 25, 2012
 *      Author: joohwile
 */

#include <iostream>

#include "MatrixCode.h"
#include "GramSchmidt.h"
#include "ConjGrad.h"
#include "Jacobi.h"
#include "OptiCode.h"
#include "SVM.h"

using namespace std;
using namespace MathCode;

void testGramSchmidt() {
	float dataA[9] = { 1, 1, -1, 1, 0, 2, 2, -2, 3 };
	SquareMatrixR<float, 3> m, n, nT, o;
	gramschmidt(m, n);
	n.transpose(nT);
	n.mult(nT, o);
	printVar(o);

	SquareMatrixR<float, 3> A, R, Q, Qt, I;
	A.fillR(dataA);
	gramschmidtQR<float, 3>(A, Q, R);
	Q.transpose(Qt);
	Q.mult(Qt, I);
	printVar(A);
	printVar(Q);
	printVar(Qt);
	printVar(R);
	printVar(I);

	gramschmidtQRx(A, Q, R);
	printVar(A);
	printVar(Q);
	printVar(Qt);
	printVar(R);
	printVar(I);

}

void testConjGrad() {
	float d[4] = { 4, 1, 1, 3 };

	Mat2 m;
	m.fillR(d);

	Vec2 b, x0, xs;
	b.set(1, 2);
	x0.set(2, 1);

	cout << "b = " << b << endl;
	cout << "x0 = " << x0 << endl;

	ConjGrad<float, 2> cg(&m, &b, &x0, &xs);
	cg.compute();

	cout << "real solution: [0.0909, 0.6364]" << endl;
	cout << "cg solution: " << xs << endl;
}

void testJacobi() {
	float dataA[9] = { 1, 1, -1, 1, 0, 2, 2, -2, 3 };

	Mat3 A;
	A.fillR(dataA);
	Vec3 b(1, 0, 0), x0(1, 0, 0), xs;

	jacobi_solve(&A, &b, &x0, &xs);
	cout << "Solution = " << xs << endl;

}

void testMult() {
	Mat3 a, b, r, c;
	float xa[9] = { 3, -1, 0, 2, 5, 1, -7, 1, 3 };
	float xb[9] = { 6, -1, 0, 0, 1, -2, 3, -8, 1 };
	a.fillR(xa);
	b.fillR(xb);
	a.mult(b, r);
	a.multC(b, c);
	cout << r << endl;
	cout << c << endl;
}

int main() {
	//testGramSchmidt();
	//testJacobi();
	testSVM();
	return 0;
}
