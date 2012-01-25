package unc;

import java.util.ArrayList;

import processing.core.PApplet;

public class Main extends PApplet {

	ArrayList pointSet;

	boolean showVoronoi = true;
	boolean showFurthestVoronoi = true;
	boolean showIntersections = false;
	boolean showConvexHull = false;
	private boolean showAnnulus = true;

	private FVoronoi fv;
	private CVoronoi cv;
	private Annulus annulus;


	public void setup() {
		pointSet = new ArrayList();

		size(500, 500);
		frameRate(30);
	}

	public void draw() {
		smooth();
		background(255);

		drawSites();
		drawVoronoi();
		drawFurthestVoronoi();
		drawConvexHull();
		drawIntersections();
		drawAnnulus();
	}

	private void drawIntersections() {
		if (showIntersections && annulus != null && annulus.is != null
				&& annulus.is.intersections != null) {
			for (int j = 0; j < annulus.is.intersections.size(); j++) {
				stroke(102,0,51);
				strokeWeight(3);
				fill(102,0,51);
				Node v = (Node) annulus.is.intersections.get(j);
				ellipse(v.coord.x, v.coord.y, 5, 5);
			}
		}
	}

	private void drawConvexHull() {
		if (showConvexHull) {
			if (fv != null && fv.inputCH != null) {
				for (int j = 0; j < fv.inputCH.size() - 1; j++) {
					stroke(90, 90, 255, 100);
					strokeWeight(3);

					QVector v1 = (QVector) fv.inputCH.get(j);
					QVector v2 = (QVector) fv.inputCH.get(j + 1);
					line(v1.x, v1.y, v2.x, v2.y);
				}
			}
		}
	}

	private void drawSites() {
		for (int j = 0; j < pointSet.size(); j++) {
			stroke(31, 120, 180);
			strokeWeight(4);
			fill(166, 206, 227);

			QVector v = (QVector) pointSet.get(j);
			ellipse((float) v.x, (float) v.y, 7, 7);
		}
	}

	private void drawVoronoi() {
		if (showVoronoi && cv != null) {
			for (int j = 0; j < cv.cvEdges.size(); j++) {
				Edge v = (Edge) cv.cvEdges.get(j);
				stroke(255, 99, 99);
				strokeWeight(1);
				line((float) v.origin.coord.x, (float) v.origin.coord.y,
						(float) v.dest.coord.x, (float) v.dest.coord.y);
			}
			for (int j = 0; j < cv.cvNodes.size(); j++) {
				Node v = (Node) cv.cvNodes.get(j);
				ellipse(v.coord.x, v.coord.y, 3, 3);
			}
		}
	}

	private void drawFurthestVoronoi() {
		if (showFurthestVoronoi) {
			if (fv == null || fv.anchors == null) {
				return;
			}

			for (int i = 0; i < fv.fvEdges.size(); i++) {
				Edge e = (Edge) fv.fvEdges.get(i);
				stroke(110,147,37);
				strokeWeight(1);
				line(e.origin.coord.x, e.origin.coord.y, e.dest.coord.x,
						e.dest.coord.y);
			}

			for (int i = 0; i < fv.fvNodes.size(); i++) {
				Node n = (Node) fv.fvNodes.get(i);
				fill(90, 200, 90);
				ellipse(n.coord.x, n.coord.y, 3, 3);
			}
		}
	}



	public void keyPressed() {

		if (keyCode == 73) {
			// i
			showIntersections = !showIntersections;
		} else if (keyCode == 86) {
			// v
			showVoronoi = !showVoronoi;
		} else if (keyCode == 70) {
			// f
			showFurthestVoronoi = !showFurthestVoronoi;
		} else if (keyCode == 67) {
			showConvexHull = !showConvexHull;
		} else if (keyCode == 65) {
			showAnnulus = !showAnnulus;
		}
	}

	public void clear() {
		pointSet.clear();
	}

	public void mousePressed() {
		if (mouseButton == LEFT) {
			QVector pIn = new QVector();
			pIn.x = mouseX;
			pIn.y = mouseY;
			pIn.ID = pointSet.size();
			pointSet.add(pIn);

			cv = new CVoronoi(pointSet, 0, 0, width, height);
			fv = new FVoronoi(pointSet, 0, 0, width, height);

			annulus = new Annulus(cv, fv, pointSet);

		} else {
			pointSet.clear();
			cv = null;
			fv = null;
			annulus = null;
		}
	}

	public void drawTwoCircle(Annulus.TwoCircle circle) {
		fill(0, 0, 0, 0);
		ellipse(circle.c.x, circle.c.y, 2 * circle.innerR, 2 * circle.innerR);
		ellipse(circle.c.x, circle.c.y, 2 * circle.outerR, 2 * circle.outerR);
	}

	public void drawAnnulus() {
		if (annulus != null && showAnnulus) {
			strokeWeight(3);
			stroke(0,26,51);
			if (annulus.best != null) {
				drawTwoCircle(annulus.best);
				fill(0, 0, 0, 255);
				text("Minimum Width = " + annulus.best.width + "    " + "Maximum Width = " + annulus.worst.width, 10, 15);
			}
			stroke(60, 60, 60, 100);
			strokeWeight(1);
			if (annulus.worst != null) {
				drawTwoCircle(annulus.worst);
			}
		}
	}

	/**
	 * 
	 */
	private static final long serialVersionUID = 3206847208968227199L;

}