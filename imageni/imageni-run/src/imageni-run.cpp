//============================================================================
// Name        : imageni-run.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "Image.h"
#include "jpegio.h"
#include "LinearTransform.h"
#include "SSDFunc.h"

using namespace std;

void visualizeSSD() {
	IntImage frame, crop;

	loadImage("", frame);
	loadImage("", crop);

	TransformParamType t;
	float v = 0;

	CSSD<int> ssd(&frame, &crop);
	ssd.position(t);
	ssd.value(v);

}

int main() {
	CImage<int> image, tII;
	image.createEmptyOf(5, 5);
	for (int i = 0; i < 25; i++) {
		image[i] = 1;
	}

	tII.createTiltedIntegralImageOf(image);
	cout << tII(0, 0) << endl;
	cout << tII(0, 1) << endl;
	cout << tII(1, 0) << endl;
	cout << tII(1, 1) << endl;
	cout << tII(3, 2) << endl;
	return 0;
}
