/*
 * testSVM.cpp
 *
 *  Created on: Jul 30, 2012
 *      Author: joohwi
 */

#include "SVM.h"
#include <iostream>

using namespace MathCode;

void testSVM() {
	typedef SVM<2> SVMType;

	typename SVMType::InputType svmInput[5];
	svmInput[0].set(1, 1);
	svmInput[1].set(2, 1);
	svmInput[2].set(2, 5);
	svmInput[3].set(-1, -1);
	svmInput[4].set(-1, -5);

	int Y[5] = { 1, 1, 1, -1, -1 };
	SVMType svm(svmInput, Y, 5);
	svm.printYXij(std::cout);

}

