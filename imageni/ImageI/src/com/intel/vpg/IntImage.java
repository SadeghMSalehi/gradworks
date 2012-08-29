package com.intel.vpg;

import java.nio.IntBuffer;

public class IntImage implements ImageBase {
	public final int _width;
	public final int _height;
	public final int[] data;

	public IntImage(int w, int h) {
		_width = w;
		_height = h;
		data = new int[w * h];
	}

	public IntImage(int[] data, int w, int h) {
		this.data = data;
		this._width = w;
		this._height = h;
	}

	public int capacity() {
		return _width * _height;
	}

	public IntBuffer asIntBuffer() {
		return IntBuffer.wrap(data);
	}

	@Override
	public int getWidth() {
		return _width;
	}

	@Override
	public int getHeight() {
		return _height;
	}

	@Override
	public int[] getIntArray() {
		return data;
	}

	@Override
	public double[] getDoubleArray() {
		return null;
	}

	public IntImage crop(Rect rect) {
		int[] data = new int[rect.w * rect.h];
		for (int j = rect.y; j < rect.y + rect.h; j++) {
			System.arraycopy(this.data, rect.x + j * _width, data, (j - rect.y)
					* rect.w, rect.w);
		}
		return new IntImage(data, rect.w, rect.h);
	}

	public IntImage createHalf() {
		int w = _width / 2;
		int h = _height / 2;
		return resize(w, h);
	}

	public IntImage resize(int w, int h) {
		IntImage resized = new IntImage(w, h);
		NativeCode.resize(data, _width, _height, resized.data, w, h, false);
		return resized;
	}

	public String toString() {
		return "image[" + _width + ", " + _height + "]";
	}

	public IntImage convolution(float[] mask) {
		IntImage pass1 = new IntImage(_width, _height);
		for (int j = 0; j < _height; j++) {
			for (int i = 0; i < _width; i++) {
				float r = 0;
				for (int k = 0; k < mask.length; k++) {
					int l = i + (k - mask.length / 2);
					if (l < 0 || l >= _width) {
						continue;
					}
					r += data[l + j * _width] * mask[k];
				}
				pass1.data[i + j * _width] = Math.round(r);
			}
		}
		IntImage pass2 = new IntImage(_width, _height);
		for (int j = 0; j < _width; j++) {
			for (int i = 0; i < _height; i++) {
				float r = 0;
				for (int k = 0; k < mask.length; k++) {
					int l = i + (k - mask.length / 2);
					if (l < 0 || l >= _height) {
						continue;
					}
					r += pass1.data[l * _width + j] * mask[k];
				}
				pass2.data[j + i * _width] = Math.round(r);
			}
		}
		return pass2;
	}

	public IntImage gaussianBlur(float sigma) {
		float[] mask = new float[7];
		float sum = 0;
		for (int i = 0; i < mask.length; i++) {
			float x = i - (mask.length / 2);
			mask[i] = (float) (Math.exp(-(x * x) / (sigma * sigma)));
			sum += mask[i];
		}
		for (int i = 0; i < mask.length; i++) {
			mask[i] /= sum;
		}

		return convolution(mask);
	}

}
