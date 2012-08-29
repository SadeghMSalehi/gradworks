package com.intel.vpg;

public class MotionParameter {
	public float tx;
	public float ty;
	public float thetaInDegree;
	public float scalePercentage;
	public Matrix forwardTransform;
	public Matrix backwardTransform;
	public String tt;
	
	public MotionParameter() {
		tx = 0;
		ty = 0;
		thetaInDegree = 0;
		scalePercentage = 100;
	}
	
	
	
	public String toString() {
		return String.format("[%f, %f, %f, %f]", tx, ty, thetaInDegree, scalePercentage);
	}

	public void copyFrom(MotionParameter x0) {
		tx = x0.tx;
		ty = x0.ty;
		thetaInDegree = x0.thetaInDegree;
		scalePercentage = x0.scalePercentage;
	}

	public float[] params() {
		return new float[] { tx, ty, thetaInDegree, scalePercentage };
	}

	public void set(float[] s) {
		tx = s[0];
		ty = s[1];
		thetaInDegree = s[2];
		scalePercentage = s[3];
	}
}
