package com.intel.vpg;

/**
 * This Face Tracker asynchronously detects a face from a given image frame.
 * When detection succeeds, this calls onDetectionSuccess() When tracking
 * succeeds, onTrackingSuccess() will be called
 * 
 * @author joohwile
 * 
 */
public class FaceTracker extends Thread {
	public static interface TrackingListener {
		public void onDetectionSuccess(FaceTracker sender, Rect faceCoord,
				IntImage detectedFace);

		public void onTrackingSuccess(FaceTracker sender,
				MotionParameter motion, IntImage trackedFace);

		public void onDetectionStart(FaceTracker faceTracker);
	}

	public enum Status {
		WAITING, DETECTING, TRACKING, DEAD
	}

	private Status _status;
	private TrackingContext _context;
	private IntImage _trackingFace;
	private TrackingListener _trackingListener;
	private IntImage _inputFrame;

	private boolean _running;
	private Rect _detectionCoord;
	private MotionParameter _currentMotion;
	private IntImage _mask;
	private IntImage _weights;

	public FaceTracker() {
		_running = true;
		_status = Status.WAITING;
		_context = new TrackingContext();
	}

	public synchronized void stopRunning() {
		this._running = false;
		notifyAll();
	}

	public synchronized void shutdownAndWait() {
		stopRunning();
		while (_status != Status.DEAD) {
			try {
				wait();
			} catch (InterruptedException e) {
			}
		}
	}
	
	public Status getStatus() {
		return _status;
	}

	public void setPixelWeights(IntImage mask) {
		_mask = mask;
	}
	
	public void setTrackingListener(TrackingListener listener) {
		_trackingListener = listener;
	}

	public synchronized void waitForDetectionReady() {
		while (_trackingFace == null) {
			try {
				wait();
			} catch (InterruptedException e) {
			}
		}
	}

	public void waitForTrackingDone() {
		while (_currentMotion == null) {
			try {
				wait();
			} catch (InterruptedException e) {
			}
		}
	}

	/**
	 * Provide new frame to track an object This frame may be ignored, if
	 * feeding occurs before tracking is done.
	 * 
	 * @param frame
	 */
	public synchronized void feedFrame(IntImage frame) {
		_inputFrame = frame;
		this.notifyAll();
	}
	
	/**
	 * Provide new frame and wait until the frame is processed.
	 * This enables synchronous tracking
	 * 
	 * @param frame
	 */
	public synchronized void feedFrameAndWait(IntImage frame) {
		if (frame == null) {
			return;
		}
		feedFrame(frame);
		while (_inputFrame != null) {
			try {
				wait();
			} catch (InterruptedException e) {
			}
		}
	}


	/**
	 * Take frame to track information
	 * 
	 * @return
	 */
	private synchronized IntImage takeFrame() {
		while (_inputFrame == null && _running) {
			try {
				_status = Status.WAITING;
				wait();
			} catch (InterruptedException e) {
			}
		}
		IntImage currentFrame = _inputFrame;
		_inputFrame = null;
		notifyAll();
		return currentFrame;
	}

	/**
	 * Asynchronous face tracking.
	 * 
	 * If new frame is not available, wait until new frame is ready If previous
	 * tracking result is not satisfying, re-run face detection (initialization)
	 * If previous tracking is good, keep tracking with new frame
	 */
	@Override
	public void run() {
		boolean keepTracking = false;
		while (_running) {
			IntImage currentFrame = takeFrame();
			if (currentFrame == null) {
				continue;
			}
			if (!keepTracking) {
				_status = Status.DETECTING;
				_trackingListener.onDetectionStart(this);
				if (beginTracking(currentFrame)) {
					keepTracking = true;
					_trackingListener.onDetectionSuccess(this, _detectionCoord,
							_trackingFace);
				}
			} else {
				_status = Status.TRACKING;
				keepTracking = _context.track(currentFrame);
				if (keepTracking) {
					setCurrentMotion(_context.getCurrentTransform());
					_trackingListener.onTrackingSuccess(FaceTracker.this,
							_currentMotion, _context.capture());
				} else {
					_currentMotion = null;
				}
			}
		}
		synchronized (this) {
			_status = Status.DEAD;
			notifyAll();
		}
	}

	private synchronized void setCurrentMotion(MotionParameter motion) {
		_currentMotion = motion;
		notifyAll();
	}

	private FaceCoords detectFaces(IntImage image, int minScale) {
		int[] faces = NativeCode.detectFaces(image, minScale);
		return new FaceCoords(faces);
	}

	/**
	 * face detection and setup for tracking
	 * 
	 * @param frame
	 * @return
	 */
	private boolean beginTracking(IntImage frame) {
		FaceCoords coords = detectFaces(frame, 2);
		if (coords == null || coords.size() == 0) {
			reset();
			return false;
		}
		
		_detectionCoord = coords.getRect(0);
		setTrackingFace(frame);
		if (_mask != null) {
			_weights = _mask.resize(_detectionCoord.w, _detectionCoord.h);
		}

		MotionParameter x0 = new MotionParameter();
		x0.tx = _detectionCoord.x + _detectionCoord.w / 2;
		x0.ty = _detectionCoord.y + _detectionCoord.h / 2;
		x0.thetaInDegree = 0;
		x0.scalePercentage = 100;
		_context.init(_trackingFace, x0, _weights);
		return true;
	}

	private synchronized void setTrackingFace(IntImage frame) {
		IntImage crop = frame.crop(_detectionCoord);
		_trackingFace = crop.gaussianBlur(3);
		notifyAll();
	}

	private void reset() {
		_detectionCoord = null;
		_trackingFace = null;
	}


}
