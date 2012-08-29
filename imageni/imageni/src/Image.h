/*
 * ImageFunction.h
 *
 *  Created on: Jun 27, 2012
 *      Author: joohwile
 */
#ifndef IMAGEFUNCTION_H_
#define IMAGEFUNCTION_H_

#include "config.h"
#include <vector>
#include <iostream>
#include <limits.h>
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include "MathCode.h"

using namespace std;
using namespace cmath;

#define __MIN(x,y)((x<y)?x:y)
#define __MAX(x,y)((x>y)?x:y)
#define __LOG(...) printf(__VA_ARGS__);

#define MAX_DBL 1e9
#define MIN_DBL 1e9

template<typename T>
class CImageTools;

template<typename T>
class CImage {
private:
	T* _ptr;
	int _width;
	int _height;
	int _bytes;

	// sometimes the memory layout is not match with _width * _height
	// for such a case it is necessary to set offset and stride explicitly
	// for most cases _offset = 0 and _stride = _width
	int _offset;
	int _stride;
	bool _deallocateMemory;

	friend class CImageTools<T> ;

	inline void nullify() {
		_width = _height = _offset = _stride = 0;
		_ptr = NULL;
	}

	void allocate() {
		_bytes = _width * _height;
		_ptr = new T[_bytes];
		if (_ptr == NULL) {
			nullify();
		} else {
			_deallocateMemory = true;
		}
	}

	inline bool allocateBuffer(int w, int h) {
		if (_ptr || !_deallocateMemory) {
			return false;
		}
		_width = w;
		_height = h;
		_offset = 0;
		_stride = _width;
		allocate();
		return true;
	}

	void deallocateBuffer() {
		if (_ptr && _deallocateMemory) {
			delete[] _ptr;
			nullify();
		}
	}

	void copyReference(CImage<T>& donor) {
		reset();
		_ptr = donor._ptr;
		_width = donor._width;
		_height = donor._height;
		_bytes = _width * _height;
		_stride = donor._stride;
		_offset = donor._offset;
	}

	inline int getHalfWidth() const {
		return _width / 2;
	}

	inline int getHalfHeight() const {
		return _height / 2;
	}

public:
	CImage() :
			_deallocateMemory(true) {
		nullify();
	}

	CImage(int w, int h) {
		_width = w;
		_height = h;
		_offset = 0;
		_stride = _width;
		allocate();
	}

	CImage(T* ptr, int w, int h) :
			_deallocateMemory(false) {
		_ptr = ptr;
		_width = w;
		_height = h;
		_offset = 0;
		_stride = _width;
		_bytes = _width * _height * sizeof(T);
	}

	/**
	 * create an empty image with same size of the given image
	 */
	template<typename U>
	CImage(CImage<U>& src) {
		allocateBuffer(src.getWidth(), src.getHeight());
	}

	virtual ~CImage() {
		deallocateBuffer();
	}

	inline int sizeOfPixel() const {
		return sizeof(T);
	}

	void setDeallocationMemory(bool b) {
		_deallocateMemory = b;
	}

	inline void zero() {
		memset(_ptr, 0, _height * _width * sizeof(T));
	}
	inline bool reset() {
		if (!_deallocateMemory) {
			return false;
		}
		deallocateBuffer();
		return true;
	}

	inline bool isEmpty() const {
		return _ptr == NULL;
	}

	inline bool isContiguous() const {
		return _offset == 0 && _stride == _width;
	}

	template<typename U>
	inline bool isSameSize(CImage<U>& ref) {
		return _width == ref.getWidth() && _height == ref.getHeight();
	}

	inline bool isInside(double x, double y) const {
		return x >= 0 && x < _width && y >= 0 && y < _height;
	}

	void wrappingOf(T* ptr, int w, int h) {
		if (_ptr != NULL) {
			return;
		}
		_ptr = ptr;
		_width = w;
		_height = h;
		_offset = 0;
		_stride = w;
		_bytes = sizeof(T) * w * h;
		_deallocateMemory = false;
	}

	template<typename U>
	bool createCopyOf(U* ptr, int w, int h) {
		if (!_deallocateMemory) {
			return false;
		}
		allocateBuffer(w, h);
		for (int i = 0; i < getNumberOfPixels(); i++) {
			_ptr[i] = static_cast<T>(ptr[i]);
		}
		return true;
	}

	bool createEmptyOf(int w, int h) {
		if (isEmpty()) {
			allocateBuffer(w, h);
			return true;
		}
		return false;
	}

	bool reuseOrCreateEmptyOf(int w, int h) {
		if (isEmpty()) {
			allocateBuffer(w, h);
			return true;
		} else if (_width == w && _height == h) {
			return true;
		}
		return false;
	}

	template<typename U>
	void createCopyOf(CImage<U>& im) {
		if (im.isContiguous()) {
			createCopyOf(im.getData(), im.getWidth(), im.getHeight());
		} else {
			if (allocateBuffer(im.getWidth(), im.getHeight())) {
				for (int j = 0; j < _height; j++) {
					T* dst = this->getScanline(j);
					U* src = im.getScanline(j);
					memcpy(dst, src, sizeof(U) * _stride);
				}
			}
		}
	}

	// nearest neighbor resize
	void createHalfOf(CImage &img, int interp = 0) {
		allocateBuffer(img.getHalfWidth(), img.getHalfHeight());
		if (interp == 0) {
			// nearest neighbor
			for (int j = 0; j < _height; j++) {
				for (int i = 0; i < _width; i++) {
					set(i, j, img(2 * i, 2 * j));
				}
			}
		} else if (interp == 1) {

		}
	}

	void createCropOf(CImage<T>& im, Vec2 xy, Vec2 wh) {
		int x = xy[0], y = xy[1];
		allocateBuffer(wh[0], wh[1]);
		for (int j = 0; j < this->getHeight(); j++) {
			T* dst = this->getScanline(j);
			const T* src = im.getScanline(j + y) + x;
			memcpy(dst, src, im.sizeOfPixel() * _stride);
		}
	}

	void createCropOf(CImage<T>& im, int coords[4]) {
		int x = coords[0], y = coords[1];
		allocateBuffer(coords[2], coords[3]);

		for (int j = 0; j < this->getHeight(); j++) {
			T* dst = this->getScanline(j);
			const T* src = im.getScanline(j + y) + x;
			memcpy(dst, src, im.sizeOfPixel() * _stride);
		}
	}

	/*
	 *
	 *
	 */
	void createResized(CImage &img, int newWidth, int newHeight,
			int method = 1) {
		if (!isEmpty()) {
			if (_width != newWidth && _height != newHeight) {
				return;
			}
		} else {
			allocateBuffer(newWidth, newHeight);
		}

		float wr = ((float) newWidth / img._width);
		float hr = ((float) newHeight / img._height);

		float ux = 1 / wr, uy = 1 / hr;

		if (method == 1) {
			for (int j = 0; j < _height; j++) {
				for (int i = 0; i < _width; i++) {
					int w = 0;
					int pc = img.bilinear(i * ux, j * uy);
					int pw = img.getAndInc((i - 1) * ux, j * uy, w);
					int pe = img.getAndInc((i + 1) * ux, j * uy, w);
					int pn = img.getAndInc(i * ux, (j - 1) * uy, w);
					int ps = img.getAndInc(i * ux, (j + 1) * uy, w);
					set(i, j,
							cmath::round(
									(pw + pe + pn + ps + 4 * pc)
											/ (float) ((w + 4))));
				}
			}
		} else if (method == 0) {
			for (int j = 0; j < _height; j++) {
				for (int i = 0; i < _width; i++) {
					int pc = img(cmath::round(i * ux), cmath::round(j * uy));
					set(i, j, pc);
				}
			}
		}
	}

	// out image must be larger than current image
	bool createPaddingOf(CImage<T>& inImg) {
		// only works with empty image
		if (_ptr) {
			return false;
		}
		_width = getNearestTwoSquares(inImg.getWidth());
		_height = getNearestTwoSquares(inImg.getHeight());
		allocateBuffer(_width, _height);
		int dw0 = (_width - inImg._width) / 2;
		int dw1 = inImg._width + dw0;
		int dh0 = (_height - inImg._height) / 2;
		int dh1 = inImg._height + dh0;
		for (int j = 0; j < _height; j++) {
			for (int i = 0; i < _width; i++) {
				if (i >= dw0 && j >= dh0 && i < dw1 && j < dh1) {
					set(i, j, inImg(i - dw0, j - dh0));
				}
			}
		}
		return true;
	}

	void giveAwayTo(CImage<T>& o) {
		o.copyReference(*this);
		nullify();
	}

	void fillIn(T value) {
		for (int i = 0; i < getNumberOfPixels(); i++) {
			_ptr[i] = value;
		}
	}

	void multiply(T coef) {
		for (int i = 0; i < _width * _height; i++) {
			_ptr[i] = _ptr[i] * coef;
		}
	}

	void inc(int x, int y) {
		set(x, y, (*this)(x, y) + 1);
	}

	inline int getNumberOfPixels() const {
		return _width * _height;
	}

	inline int getStride() const {
		return _stride;
	}

	inline int getOffset() const {
		return _offset;
	}

	inline int getWidth() const {
		return _width;
	}

	inline int getHeight() const {
		return _height;
	}

	inline int getNearestTwoSquares(int x) {
		return (int) (powf(2.0f, floorf(log10f(x) / log10f(2)) + 1));
	}

	inline void set(int x, int y, T v) {
		_ptr[x + y * _stride + _offset] = v;
	}

	inline T& operator()(int x, int y) const {
		return _ptr[x + y * _stride + _offset];
	}

	/**
	 * nearest-neighbor
	 */
	inline T& operator()(double jx, double jy) const {
		int x = (int) (jx + .5);
		int y = (int) (jy + .5);
		return _ptr[x + y * _stride + _offset];
	}

	inline T bilinear(double jx, double jy) const {
		int x0 = (int) jx;
		int y0 = (int) jy;
		double dx = jx - x0, dx1 = 1 - dx;
		double dy = jy - y0, dy1 = 1 - dy;
		int base = x0 + _stride * y0 + _offset;
		int i00 = _ptr[base];
		int i10 = _ptr[base + 1];
		int i01 = _ptr[base + _stride];
		int i11 = _ptr[base + _stride + 1];
		return dx1 * dy1 * i00 + dx * dy1 * i10 + dx1 * dy * i01 + dx * dy * i11;
	}

	inline float bilinearFloat(double jx, double jy) const {
		int x0 = (int) jx;
		int y0 = (int) jy;
		float dx = jx - x0, dx1 = 1 - dx;
		float dy = jy - y0, dy1 = 1 - dy;
		int base = x0 + _stride * y0 + _offset;
		int i00 = _ptr[base];
		int i10 = _ptr[base + 1];
		int i01 = _ptr[base + _stride];
		int i11 = _ptr[base + _stride + 1];
		return dx1 * dy1 * i00 + dx * dy1 * i10 + dx1 * dy * i01 + dx * dy * i11;
	}

	inline float safe_bilinear(double jx, double jy) const {
		int x0 = (int) jx;
		int y0 = (int) jy;
		double dx = jx - x0, dx1 = 1 - dx;
		double dy = jy - y0, dy1 = 1 - dy;
		int base = x0 + _stride * y0 + _offset;
		if (base >= 0 && base + _stride + 1 < _width * _height) {
			int i00 = _ptr[base];
			int i10 = _ptr[base + 1];
			int i01 = _ptr[base + _stride];
			int i11 = _ptr[base + _stride + 1];
			return dx1 * dy1 * i00 + dx * dy1 * i10 + dx1 * dy * i01
					+ dx * dy * i11;
		} else {
			return 0;
		}
	}

	inline T getAndInc(float x, float y, int& w) const {
		if (x >= 0 && x < _width && y >= 0 && y < _height) {
			w++;
			return cmath::round((*this).bilinearFloat(x, y));
		}
		return 0;
	}

	inline T& operator[](int i) const {
		return _ptr[i + _offset];
	}

	T getSumRect(int x, int y, int w, int h) const {
		typedef T U;
		const int base = y * _stride + _offset, tl = base + x, tr = tl + w, bl =
				tl + h * _width, br = bl + w;
		return (U) (_ptr[tl] - _ptr[tr] - _ptr[bl] + _ptr[br]);
	}

	inline T minValue() const {
		double min = MIN_DBL;
		if (isContiguous()) {
			for (int i = 0; i < _width * _height; i++) {
				min = __MIN(min, _ptr[i]);
			}
		} else {
			for (int j = 0; j < _height; j++) {
				T* pos = _ptr + _offset + (j * _stride);
				for (int i = 0; i < _width; i++) {
					min = __MIN(min, pos[i]);
				}
			}

		}
		return (T) min;
	}

	inline T maxValue() const {
		double max = -MAX_DBL;
		for (int j = 0; j < _height; j++) {
			T* pos = this->getScanline(j);
			for (int i = 0; i < _width; i++) {
				max = __MAX(max, pos[i]);
			}
		}

		return (T) max;
	}

	T* getData() const {
		return _ptr + _offset;
	}

	T* getScanline(int y) const {
		return &_ptr[_offset + y * _stride];
	}

	template<typename U, typename V>
	inline void dot(CImage<U>& uImage, CImage<V>& vImage) {
		if (!uImage.isEmpty() && !vImage.isEmpty()) {
			int N = _width * _height;
			for (int i = 0; i < N; i++) {
				vImage._ptr[i] = _ptr[i] * uImage._ptr[i];
			}
		}
	}

	template<typename U>
	void createRescaled(CImage<U> &inImg) {
		if (isEmpty()) {
			if (!allocateBuffer(inImg.getWidth(), inImg.getHeight())) {
				return;
			}
		}
		if (_width == inImg.getWidth() && _height == inImg.getHeight()) {
			double maxT = inImg.maxValue();
			double minT = inImg.minValue();
			for (int j = 0; j < _height; j++) {
				const U* src = inImg.getScanline(j);
				T* dst = this->getScanline(j);
				for (int i = 0; i < _width; i++) {
					dst[i] = 255 * ((src[i] - minT) / maxT);
				}
			}
			return;
		}
		throw "Destination not available!!";
	}

	template<typename ReturnType, typename MaskType>
	inline ReturnType computeFilterResponse(int ii, int jj,
			CImage<MaskType>& mask) {
		if (isContiguous()) {
			int radX = mask.getWidth() / 2;
			int radY = mask.getHeight() / 2;
			int w0 = __MAX(ii-radX,0);
			int w1 = __MIN(ii+radX+1,_width);
			int h0 = __MAX(jj-radY,0);
			int h1 = __MIN(jj+radY+1,_height);
			ReturnType response = 0;
			for (int j = h0, cnt = 0; j < h1; j++) {
				for (int i = w0; i < w1; i++) {
					response = response
							+ (_ptr[i + j * _width] * mask.getData()[cnt++]);
				}
			}
			return response;
		} else {
			return -1;
		}
	}

	template<typename U>
	void createSquaredOf(CImage<U>& inImg) {
		if (!allocateBuffer(inImg.getWidth(), inImg.getHeight())) {
			return;
		}
		for (int j = 0; j < _height; j++) {
			U* src = inImg.getScanline(j);
			T* dst = this->getScanline(j);
			for (int i = 0; i < _width; i++) {
				dst[i] = src[i] * src[i];
			}
		}
	}

	template<typename U>
	inline bool computeGaussian(CImage<U> &out, int w, int sigma) {
		if (!out.createEmpty(*this)) {
			return false;
		}
		if (!this->isContiguous()) {
			return false;
		}
		// w must be odd number
		float* maskData = new float[w];
		int r = w / 2;
		for (int x = -r; x <= r; x++) {
			maskData[x + r] = expf(-(x * x) / 2 / sigma / sigma);
		}
		float sum = 0;
		for (int x = 0; x < w; x++) {
			sum += maskData[x];
		}
		for (int x = 0; x < w; x++) {
			maskData[x] /= sum;
		}
		CImage<U> gX;
		if (gX.createEmpty(*this)) {
			computeConvolution(gX, maskData, w, 1);
			gX.computeConvolution(out, maskData, 1, w);
		}
		delete[] maskData;
	}

	template<typename U, typename S>
	bool computeConvolution(CImage<U>& out, CImage<S>& mask) {
		if (!this->isContiguous()) {
			return false;
		}
		U* outPtr = out.getData();
		if (outPtr == 0) {
			return false;
		}
		for (int j = 0; j < _height; j++) {
			for (int i = 0; i < _width; i++) {
				// FIXME : boundary condition
				outPtr[i + _width * j] = computeFilterResponse<U, S>(i, j,
						mask);
			}
		}
		return true;
	}

	template<typename U, typename S>
	bool computeConvolution(CImage<U>& out, S* mask, int mx, int my) {
		if (!out.createEmpty(*this)) {
			return false;
		}
		if (this->isContiguous()) {
			return false;
		}
		int rx = mx / 2;
		int ry = my / 2;
		for (int j = ry; j < _height - ry; j++) {
			for (int i = rx; i < _width - rx; i++) {
				U sum = 0;
				int cnt = 0;
				for (int y = -ry; y <= ry; y++) {
					for (int x = -rx; x <= rx; x++) {
						sum += (mask[cnt++] * (*this)(i + x, j + y));
					}
				}
				out.set(i, j, sum);
			}
		}
	}

	bool computeSobelX(CImage<T>& out) {
		T maskData[9] = { -1, 0, 1, -2, 0, 2, -1, 0, 1 };
		return computeConvolution(out, maskData, 3, 3);
	}

	bool computeSobelY(CImage<T>& out) {
		T maskData[9] = { -1, -2, -1, 0, 0, 0, 1, 2, 1 };
		return computeConvolution(out, maskData, 3, 3);
	}

	void write(ostream& os) {
		for (int j = 0; j < _height; j++) {
			os << (*this)(0, j);
			for (int i = 1; i < _width; i++) {
				os << " " << (*this)(i, j);
			}
			os << endl;
		}
	}
	/*
	 bool computeSecondOrderDerivative(CImage& out, int dir) {
	 if (out.isEmpty()) {
	 out.createEmpty(*this);
	 }
	 if (dir == 0) {
	 int maskData[9] = { 1, -2, 1, 2, -4, 2, 1, -2, 1 };
	 computeConvolution(out, maskData, 3, 3);
	 } else if (dir == 1) {
	 int maskData[9] = { 1, 2, 1, -2, -4, -2, 1, 2, 1 };
	 computeConvolution(out, maskData, 3, 3);
	 } else if (dir == 2) {
	 CImage Ix;
	 Ix.createEmpty(*this);
	 computeSobelX(Ix);
	 Ix.computeSobelY(out);
	 }
	 return true;
	 }
	 */

	template<typename U>
	void createIntegralImageOf(CImage<U> &inImg) {
		if (!allocateBuffer(inImg.getWidth(), inImg.getHeight())) {
			return;
		}
		for (int j = 0; j < _height; j++) {
			T colSum = 0;
			for (int i = 0; i < _width; i++) {
				colSum += inImg(i, j);
				if (j == 0) {
					set(i, j, colSum);
				} else {
					set(i, j, (*this)(i, j - 1) + colSum);
				}
			}
		}
	}

	template<typename U>
	void createTiltedIntegralImageOf(CImage<U>& inImg) {
		if (!allocateBuffer(inImg.getWidth(), inImg.getHeight())) {
			return;
		}

		// first pass: computes diagonal top left
		zero();
		T* out = getScanline(0);
		U* in = inImg.getScanline(0);
		out[0] = in[0];
		for (int i = 1; i < _width; i++) {
			out[i] = out[i - 1] + in[i];
		}
		for (int j = 1; j < _height; j++) {
			out = getScanline(j);
			T* outUp = getScanline(j - 1);
			in = inImg.getScanline(j);
			U* inUp = inImg.getScanline(j - 1);
			out[0] = in[0];
			out[1] = out[0] + outUp[0] + in[1];
			for (int i = 2; i < _width; i++) {
				out[i] = out[i - 1] + outUp[i - 1] + in[i] - outUp[i - 2];
			}
		}

		// second pass: computes diagonal from bottom right
		for (int j = _height - 2; j >= 0; j--) {
			out = getScanline(j);
			T* outDown = getScanline(j + 1);
			in = inImg.getScanline(j);
			for (int i = _width - 1; i > 1; i--) {
				out[i] = out[i] + outDown[i - 1] - out[i - 2];
			}
		}
	}

	/*
	 void computeSpatialMaximum(CImage<T> &out, int windowSize) {
	 if (!allocateBuffer(inImg.getWidth(), inImg.getHeight())) {
	 return;
	 }
	 int radius = windowSize / 2;

	 // ignore boundary
	 for (int j = radius; j < _height - radius; j++) {
	 for (int i = radius; i < _width - radius; i++) {
	 T maxResp = -INT_MAX;
	 for (int y = -radius; y <= radius; y++) {
	 for (int x = -radius; x <= radius; x++) {
	 maxResp =
	 maxResp < (*this)(i + x, j + y) ?
	 (*this)(i + x, j + y) : maxResp;
	 }
	 }
	 if (maxResp > (*this)(i, j)) {
	 out.set(i, j, 0);
	 } else {
	 out.set(i, j, maxResp);
	 }
	 }
	 }
	 }
	 */

};

#ifdef JNI
#include "jni.h"
typedef CImage<jint> IntImage;
typedef CImage<jfloat> FloatImage;
typedef CImage<jdouble> DoubleImage;
#else
typedef CImage<int> IntImage;
typedef CImage<float> FloatImage;
typedef CImage<double> DoubleImage;
#endif

#endif /* IMAGEFUNCTION_H_ */

