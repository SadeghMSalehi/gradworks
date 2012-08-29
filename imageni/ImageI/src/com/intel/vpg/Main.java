package com.intel.vpg;

import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import com.intel.vpg.FaceTracker.TrackingListener;

public class Main {
	public static final String DATA_DIR = "../";
//	public static final String INPUT_SEQUENCES = DATA_DIR
//			+ "data/myface/hard/myface-%03d.jpg";
	public static final String INPUT_SEQUENCES = "c:/Users/joohwile/Desktop/tracking2/face%04d.jpg";
	public static final String OUTPUT_FOLDER = "c:/Users/joohwile/Desktop/tracking2/output";
	public static final int NUM_FRAMES = 1000;

	public static void main(String[] args) {
		// testBlur();
		//singleTracking();
		multipleTracking();
		// emmaWatson();
	}

	private static void detectionTest() {
		IntImage frame = readImage("");
		NativeCode.detectFaces(frame, 12);
	}
	private static void testBlur() {
		IntImage frame0 = readImage(DATA_DIR + "/data/myface/hard/myface-001.jpg");
		IntImage blur = frame0; // frame0.gaussianBlur(3);
		writeJPEG("blur.jpg", blur);
	}

	private static void multipleTracking() {
		int iter = 0;
		FaceTracker tracker = new FaceTracker();
		tracker.setTrackingListener(new TrackingListener() {
			private int _counter = 0;

			@Override
			public void onTrackingSuccess(FaceTracker sender,
					MotionParameter motion, IntImage trackedFace) {
				writeJPEG(String.format(OUTPUT_FOLDER + "/__trackingResult_%03d.jpg",
						++_counter), trackedFace);
			}

			@Override
			public void onDetectionSuccess(FaceTracker sender, Rect faceCoord,
					IntImage detectedFace) {
				System.out.println("Detection coordinate: " + faceCoord);
				writeJPEG(String.format(OUTPUT_FOLDER + "/__detectionResult_%03d.jpg",
						++_counter), detectedFace);
			}

			@Override
			public void onDetectionStart(FaceTracker faceTracker) {
				System.out.println("Detection has started... wait for a while");
			}
		});

		tracker.start();
		for (iter = 1; iter < NUM_FRAMES; iter++) {
			String f = String.format(INPUT_SEQUENCES, iter);
			// System.out.println("Processing " + f);
			IntImage frame = readImage(f);
			if (frame == null) {
				break;
			}
			tracker.feedFrameAndWait(frame);
		}
		tracker.shutdownAndWait();
	}

	public static void singleTracking() {
		FaceTracker tracker = new FaceTracker();
		IntImage frame0 = readImage(DATA_DIR + "/data/myface/myface-001.jpg");
		IntImage frame1 = readImage(DATA_DIR + "/data/myface/myface-027.jpg");
		IntImage mask = readImage("mask.jpg");

		// tracker.setPixelWeights(mask);
		tracker.setTrackingListener(new TrackingListener() {
			@Override
			public void onTrackingSuccess(FaceTracker sender,
					MotionParameter motion, IntImage trackedFace) {
				writeJPEG(DATA_DIR
						+ "data/myface/tracking/single/__trackingResult.jpg",
						trackedFace);
			}

			@Override
			public void onDetectionSuccess(FaceTracker sender, Rect faceCoord,
					IntImage detectedFace) {
				System.out.println("Detection coordinate: " + faceCoord);
				writeJPEG(DATA_DIR
						+ "data/myface/tracking/single/__detectionResult.jpg",
						detectedFace);
			}

			@Override
			public void onDetectionStart(FaceTracker faceTracker) {
				System.out.println("Detection has started... wait for a while");
			}
		});

		tracker.start();
		tracker.feedFrameAndWait(frame0);
		tracker.feedFrameAndWait(frame1);
		tracker.shutdownAndWait();
	}

	public static void emmaWatson() {
		FaceTracker tracker = new FaceTracker();
		IntImage frame0 = readImage(DATA_DIR + "/data/emma.jpg");
		IntImage frame1 = readImage(DATA_DIR + "/data/emma_10deg.jpg");

		tracker.setTrackingListener(new TrackingListener() {
			@Override
			public void onTrackingSuccess(FaceTracker sender,
					MotionParameter motion, IntImage trackedFace) {
				System.out.println("Motion parameter: " + motion);
				writeJPEG(DATA_DIR + "/data/__trackingResult.jpg", trackedFace);
			}

			@Override
			public void onDetectionSuccess(FaceTracker sender, Rect faceCoord,
					IntImage detectedFace) {
				System.out.println("Detection coordinate: " + faceCoord);
				writeJPEG(DATA_DIR + "/data/__detectionResult.jpg", detectedFace);
			}

			@Override
			public void onDetectionStart(FaceTracker faceTracker) {
				System.out.println("Detection has started... wait for a while");
			}
		});

		tracker.start();
		tracker.feedFrameAndWait(frame0);
		tracker.feedFrameAndWait(frame1);
		tracker.shutdownAndWait();
	}
	
	/**
	 * Read image from file if no file exists, return null otherwise return
	 * gray-color converted IntImage
	 * 
	 * @param f
	 * @return
	 */
	private static IntImage readImage(String f) {
		try {
			BufferedImage rawImg = ImageIO.read(new File(f));
			IntImage image = convert(rawImg);
			IntImage halfImage = image.createHalf();
			// System.out.println("reading " + f);
			return image;
		} catch (IOException e) {
			return null;
		}
	}

	public static IntImage convert(BufferedImage image) {
		final int w = image.getWidth();
		final int h = image.getHeight();
		int[] rgbArray = new int[w * h];
		image.getRGB(0, 0, w, h, rgbArray, 0, w);
		IntImage intImg = new IntImage(image.getWidth(), image.getHeight());
		for (int i = 0; i < rgbArray.length; i++) {
			int r = (rgbArray[i] >> 16) & 0xff;
			int g = (rgbArray[i] >> 8) & 0xff;
			int b = (rgbArray[i] >> 0) & 0xff;
			int gray = (int) Math.round((r + g + b) / 3.0);
			intImg.data[i] = gray;
		}
		return intImg;
	}

	public static void writeJPEG(String file, IntImage image) {
		BufferedImage im = new BufferedImage(image.getWidth(),
				image.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
		WritableRaster raster = im.getRaster();
		raster.setPixels(0, 0, image.getWidth(), image.getHeight(), image.data);
		try {
			ImageIO.write(im, "jpg", new File(file));
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
