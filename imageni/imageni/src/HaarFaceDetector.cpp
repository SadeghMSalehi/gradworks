/*
 * HaarFaceDetector.cpp
 *
 *  Created on: Jun 29, 2012
 *      Author: joohwile
 */

#include "config.h"
#include "HaarFaceDetector.h"
#include <iostream>
#include <cmath>

using namespace std;

void loadTrainingSet(HaarCascade &trainingSet) {
#include "haar_faces.inc"
}

//template<typename T>
//inline void print(std::vector<T>& vec) {
//	for (size_t i = 0; i < vec.size(); i++) {
//		cout << vec[i] << " ";
//	}
//	cout << endl;
//}
//
//static inline int count1(std::vector<int>& vec) {
//	int cnt = 0;
//	for (size_t i = 0; i < vec.size(); i++) {
//		if (vec[i] == 1) {
//			cnt++;
//		}
//	}
//	return cnt;
//}

template <typename I, typename R>
void CHaarFaceDetector<I,R>::printTrainedInformation() {
	cout << "Number of stages: " << _trainingSet.stages.size() << endl;
	for (size_t i = 0; i < _trainingSet.stages.size(); i++) {
		cout << "Number of features: " << _trainingSet.stages[i].features.size()
				<< endl;
		for (size_t j = 0; j < _trainingSet.stages[i].features.size(); j++) {
			cout << _trainingSet.stages[i].features[j].numRects << " ";
			cout << "[" << _trainingSet.stages[i].features[j].left_val << ", "
					<< _trainingSet.stages[i].features[j].right_val << ", "
					<< _trainingSet.stages[i].features[j].threshold << "]"
					<< endl;
		}

	}
}

