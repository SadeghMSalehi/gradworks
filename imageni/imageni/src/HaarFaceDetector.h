/*
 * HaarFaceDetector.h
 *
 *  Created on: Jun 29, 2012
 *      Author: joohwile
 */

#ifndef HAARFACEDETECTOR_H_
#define HAARFACEDETECTOR_H_

#include "config.h"
#include "HaarFeature.h"
#include "Image.h"
#include <vector>

void loadTrainingSet(HaarCascade &trainingSet);

template<typename I, typename R>
class CHaarFaceDetector {
private:
	HaarCascade _trainingSet;
	haar_point_array_t _detectionResult;
	int _minScaleFactor;
public:
	CHaarFaceDetector() {
		// TODO Auto-generated constructor stub
		DEBUG_MODE = 0;
		loadTrainingSet(_trainingSet);
		_minScaleFactor = 2;
	}

	virtual ~CHaarFaceDetector() {
	}

	void setMinScaleFactor(int factor) {
		_minScaleFactor = factor;
	}

	void detectAtScales(CImage<I>& orgImg, CImage<I>& roiImg, R &detectedScale,
			bool earlyReturn = true) {
		CImage<I> II, I2, II2;
		II.createIntegralImageOf(roiImg);
		I2.createSquaredOf(roiImg);
		II2.createIntegralImageOf(I2);

		//LOGD("II computed.\n");
		const R scaleFactor = 1.2;
		const R scaleUpdate = 1 / scaleFactor;
		R scaleWidth = (float) orgImg.getWidth()
				/ (float) _trainingSet.featureWidth;
		R scaleHeight = (float) orgImg.getHeight()
				/ (float) _trainingSet.featureHeight;

		//scaleWidth /= 2;
		//scaleHeight /= 2;

		R maxScale = __MIN(scaleWidth, scaleHeight);
		R minScale = maxScale / _minScaleFactor;
		I maxIter = ceil(log(maxScale) / log(scaleFactor));
		R scale = maxScale;
		_detectionResult.clear();
		std::vector<I> x, y;

		//LOGD("Loop begins: maxIter = %d, minScale = %f\n", maxIter, minScale);
		for (int i = 0; i < maxIter && scale >= minScale; i++) {
			int fw = floor(_trainingSet.featureWidth * scale);
			int fh = floor(_trainingSet.featureHeight * scale);
			LOGD("Searching feature size = %d\n", fw);

			int step = floor(scale > 2 ? scale : 2);
			for (int j = 0; j < II.getWidth() - fw; j += step) {
				x.push_back(j);
			}
			for (int j = 0; j < II.getHeight() - fh; j += step) {
				y.push_back(j);
			}
			if (x.size() == 0 || y.size() == 0) {
				continue;
			}
			detectAtSingleScale(x, y, scale, II, fw, fh, II2);
			if (_detectionResult.size() > 3) {
				detectedScale = scale;
				if (earlyReturn) {
					break;
				}
			}
			scale = scale * scaleUpdate;
		}

	}
	haar_point_array_t getResult() {
		return _detectionResult;
	}
	int DEBUG_MODE;

private:
	void detectAtSingleScale(std::vector<I>& xx, std::vector<I>& yy, R scale,
			CImage<I>& II, I w, I h, CImage<I>& II2) {
		R inverseArea = 1.0f / (w * h);

		// sum of haar feature response at each detection location
		for (size_t iY = 0; iY < yy.size(); iY++) {
			for (size_t iX = 0; iX < xx.size(); iX++) {
				// can be parallel
				bool isDetected = true;
				I x = xx[iX];
				I y = yy[iY];
				R sum2mean = II2.getSumRect(x, y, w, h) * inverseArea;
				R mean = II.getSumRect(x, y, w, h) * inverseArea;
				R stdev = sqrtf(sum2mean - mean * mean);
				R stageSum = 0;
				for (size_t iStage = 0; iStage < _trainingSet.stages.size();
						iStage++) {
					stageSum = 0;
					HaarStage& stage = _trainingSet.stages[iStage];
					for (size_t iTree = 0; iTree < stage.features.size();
							iTree++) {
						haar_feature_t& feature = stage.features[iTree];
						float treeSum = detectTree(feature, scale, x, y, II,
								stdev, inverseArea);
						stageSum += treeSum;
					}
					if (stageSum < stage.stage_threshold) {
						isDetected = false;
						break;
					}
				}
				if (isDetected) {
					haar_point_t featurePoint;
					featurePoint.x = xx[iX];
					featurePoint.y = yy[iY];
					featurePoint.w = w;
					featurePoint.h = h;
					featurePoint.scale = scale;
					_detectionResult.push_back(featurePoint);
				}
			}
		}
	}

	R detectTree(haar_feature_t &feature, R scale, I x, I y, CImage<I>& II,
			R stdv, R inverseArea) {
		R rectSum = 0;
		for (int i = 0; i < feature.numRects; i++) {
			I rx = floor(feature.rects[i].x * scale + x);
			I ry = floor(feature.rects[i].y * scale + y);
			I rw = floor(feature.rects[i].w * scale);
			I rh = floor(feature.rects[i].h * scale);
			R weight = feature.rects[i].weight;
			R r_sum = II.getSumRect(rx, ry, rw, rh) * weight;
			rectSum += (r_sum);
		}
		rectSum = rectSum * inverseArea;
		if (rectSum >= (feature.threshold * stdv)) {
			return feature.right_val;
		}
		return feature.left_val;
	}

	void printTrainedInformation();
};

#endif /* HAARFACEDETECTOR_H_ */
