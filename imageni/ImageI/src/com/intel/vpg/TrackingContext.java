package com.intel.vpg;

import java.util.ArrayList;

public class TrackingContext {
	private IntImage _target;
	private double[] _gaussNewtonPrior;

	private MotionParameter _x0;
	private MotionParameter _xs;

	private IntImage _frame;
	private boolean _keepTracking = false;

	private ArrayList<Double> _costReport;
	private IntImage _weights;

	public boolean isTrackingReady() {
		// TODO Auto-generated method stub
		return false;
	}

	public void init(IntImage target, MotionParameter x0, IntImage weights) {
		this._xs = new MotionParameter();
		this._x0 = new MotionParameter();
		this._target = target;
		this._xs.copyFrom(x0);
		this._gaussNewtonPrior = new double[4 * target.data.length];
		this._weights = weights;
		nativeInit(target.data, target.getWidth(), target.getHeight(),
				_gaussNewtonPrior);
		if (_weights != null) {
			int stride = _weights._height * _weights._width;
			if (4 * stride == _gaussNewtonPrior.length) {
				for (int i = 0; i < stride; i++) {
					if (_weights.data[i] < 100) {
						_gaussNewtonPrior[i] = 0;
						_gaussNewtonPrior[i + stride] = 0;
						_gaussNewtonPrior[i + stride * 2] = 0;
						_gaussNewtonPrior[i + stride * 3] = 0;
					}
				}
			}
		}
	}

	public boolean track(IntImage frame) {
		_x0.copyFrom(_xs);
		_costReport = new ArrayList<Double>();
		this._frame = frame;
		float[] solution = new float[4];
		Clock clk = new Clock();
		double cost = nativeTrack(_gaussNewtonPrior, this._frame.data, frame.getWidth(),
				frame.getHeight(), _target.data, _target.getWidth(),
				_target.getHeight(), _x0.params(), 300, solution);
		if (true) {
			_keepTracking = true;
			_xs.set(solution);
			System.out.println(clk.tock() + " ms");
			return true;
		} else {
			_keepTracking = false;
			return false;
		}
	}

	public IntImage capture() {
		int w = _target.getWidth();
		int h = _target.getHeight();
		IntImage trackedTarget = new IntImage(_target.getWidth(),
				_target.getHeight());
		nativeResample(_frame.data, _frame.getWidth(), _frame.getHeight(),
				_xs.params(), trackedTarget.data, w, h);
		return trackedTarget;
	}

	public MotionParameter getCurrentTransform() {
		return _xs;
	}

	private void reportCost(double cost) {
		_costReport.add(cost);
	}

	private String getCostAsString() {
		if (_costReport == null || _costReport.size() == 0) {
			return "[]";
		}
		StringBuffer buf = new StringBuffer();
		buf.append("[").append(_costReport.get(0));
		for (int i = 0; i < _costReport.size(); i++) {
			buf.append(", ").append(_costReport.get(i));
		}
		buf.append("]");
		return buf.toString();
	}

	public native void nativeInit(int[] target, int tw, int th,
			double[] outPrior);

	public native double nativeTrack(double[] prior, int[] frame, int w,
			int h, int[] target, int tw, int th, float[] params, int maxIter,
			float[] paramsOut);

	public native boolean nativeResample(int[] frame, int w, int h,
			float[] params, int[] outTarget, int ow, int oh);

}
