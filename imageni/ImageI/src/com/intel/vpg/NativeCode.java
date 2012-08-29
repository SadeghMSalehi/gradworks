package com.intel.vpg;

public class NativeCode {
	static {
		System.loadLibrary("imageni");
	}

	public static native int[] detectFaces(IntImage image, int minScale);

	public native static void resize(int[] data, int _width, int _height, int[] data2,
			int w, int h, boolean earlyReturn);
}
