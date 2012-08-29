/*
 * ImageFeatures.h
 *
 *  Created on: Jul 3, 2012
 *      Author: joohwile
 */

#ifndef IMAGEFEATURES_H_
#define IMAGEFEATURES_H_

#include <vector>

#define BUFSIZ 1024
#define MAX_FEATURE_NUM 64
#define __LOG(...) printf(__VA_ARGS__);

class CVect;
class CFeature;
enum FeatureType {
	SIFT, MSER
};

typedef std::vector<CVect> CVectArray;
typedef std::vector<CFeature&> CFeatureArray;

class CVect {
public:
	float x, y, z;

	CVect() {
		this->x = 0;
		this->y = 0;
		this->z = 0;
	}
	CVect(float x, float y, float z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}

	CVect& operator=(const CVect& o) {
		this->x = o.x;
		this->y = o.y;
		this->z = o.z;
		return *this;
	}

	CVect& operator+(const CVect& o) {
		x += o.x;
		y += o.y;
		z += o.z;
		return *this;
	}

	CVect& operator-(const CVect& o) {
		x -= o.x;
		y -= o.y;
		z -= o.z;
		return *this;
	}

	CVect& operator*(const CVect& o) {
		x *= o.x;
		y *= o.y;
		z *= o.z;
		return *this;
	}
};

class CFeature {
private:
	int _size;
	float _features[MAX_FEATURE_NUM];
	CVect _point;
	float _angle;
	float _scale;
	FeatureType _featureType;

public:
	CFeature(int nSize, float x, float y, float angle, float scale,
			CVect point = CVect()) {
		_size = nSize;
		_angle = angle;
		_scale = scale;
		_point = point;
		if (nSize > MAX_FEATURE_NUM) {
			__LOG("Max size exceeded!");
		}
		memset(_features, 0, sizeof(_features));
	}

	CFeature(const float* features, int nSize, float x, float y, float z,
			float angle, float scale, CVect point = CVect()) {
		_size = nSize;
		_angle = angle;
		_scale = scale;
		_point = point;
		if (nSize > BUFSIZ) {
			__LOG("Max size exceeded!");
		}
		memcpy(_features, features, sizeof(_features));
	}

	virtual ~CFeature() {
	}

	int getSize() {
		return _size;
	}

	int getAngle() {
		return _angle;
	}

	int getScale() {
		return _scale;
	}

	float* getFeatures() {
		return _features;
	}

	virtual FeatureType getType() {
		return _featureType;
	}
};


template<class ImageType>
class CHarrisSIFTComputer {
private:
	float _threshold;
	int _nLayers;
	int _nMaxInterestPoints;
	float _minDistance;
	int _nLevels;
	ImageType* _pImage;
	CVectArray _interestPoints;

public:
	CHarrisSIFTComputer(float threshold = 0.01f, int nLayers = 3,
			int nMaxInterestPoints = 512) {
	}

	int ComputeFeatures(ImageType* pImage);

	int getNumberOfLevels() {
		return _nLevels;
	}
	float getMinDistance() {
		return _minDistance;
	}
	float getThreshold() {
		return _threshold;
	}
	int getMaxInterestPoints() {
		return _nMaxInterestPoints;
	}
};

#endif /* IMAGEFEATURES_H_ */
