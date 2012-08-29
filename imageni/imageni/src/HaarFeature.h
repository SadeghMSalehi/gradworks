/*
 * HaarFeature.h
 *
 *  Created on: Jun 29, 2012
 *      Author: joohwile
 */

#ifndef HAARFEATURE_H_
#define HAARFEATURE_H_

#include <vector>
#include "jni.h"

typedef struct haar_point {
	jint x;
	jint y;
	jint w;
	jint h;
	float scale;
} haar_point_t;
typedef std::vector<haar_point_t> haar_point_array_t;

typedef struct haar_rect {
	jint x,y,w,h;
	float weight;
} haar_rect_t;

typedef struct haar_feature {
	jint numRects;
	haar_rect_t rects[4];
	float threshold;
	float left_val;
	float right_val;
} haar_feature_t;

class HaarStage {
public:
	std::vector<haar_feature_t> features;
	float stage_threshold;

	void clear() {
		stage_threshold = 0;
		features.clear();
	}
};

class HaarCascade {
public:
	jint featureWidth;
	jint featureHeight;
	std::vector<HaarStage> stages;
};


#endif /* HAARFEATURE_H_ */
