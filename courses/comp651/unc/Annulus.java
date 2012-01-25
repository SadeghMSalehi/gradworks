package unc;
import java.util.ArrayList;
import java.util.Iterator;

public class Annulus {
	ArrayList case1, case2, case3;

	public class TwoCircle {
		public QVector c;
		public float innerR, outerR;
		public float width;

		public TwoCircle(QVector c, float innerR, float outerR) {
			this.c = c;
			this.innerR = innerR;
			this.outerR = outerR;
			this.width = outerR - innerR;
		}
		
		public String toString() {
			return c.toString() + " width = " + width;
		}
	}

	public TwoCircle best1, best2, best3, best;
	public TwoCircle worst1, worst2, worst3, worst;
	public Intersector is;

	public Annulus(CVoronoi cv, FVoronoi fv, ArrayList pointSet) {
		case1 = detectCase1(cv, pointSet);
		case2 = detectCase2(fv, pointSet);
		case3 = detectCase3(cv, fv, pointSet);

		if (best1 != null && best2 != null && best3 != null) {
			if (best1.width <= best2.width && best1.width <= best3.width) {
				best = best1;
			} else if (best2.width <= best3.width && best2.width <= best1.width) {
				best = best2;
			} else if (best3.width <= best1.width && best3.width <= best2.width) {
				best = best3;
			}

			if (worst1.width >= worst2.width && worst1.width >= worst3.width) {
				worst = worst1;
			} else if (worst2.width >= worst3.width
					&& worst2.width >= worst1.width) {
				worst = worst2;
			} else if (worst3.width >= worst1.width
					&& worst3.width >= worst2.width) {
				worst = worst3;
			}
		}
	}
	
	private float findMin(float[] d) {
		float min = -1;
		for (int i = 0; i < d.length; i++) {
			if (d[i] <= min || min < 0) {
				min = d[i];
			}
		}
		return min;
	}
	
	private float findMax(float[] d) {
		float max = -1;
		for (int i = 0; i < d.length; i++) {
			if (d[i] >= max || max < 0) {
				max = d[i];
			}
		}
		return max;
	}

	private ArrayList detectCase3(CVoronoi cv, FVoronoi fv, ArrayList pointSet) {
		if (fv == null || fv.fvEdges == null || cv == null
				|| cv.cvEdges == null) {
			return null;
		}

		ArrayList cand = new ArrayList();
		double minWidth = -1, maxWidth = -1;

		is = new Intersector(cv.cvEdges, fv.fvEdges);
		for (int i = 0; i < is.intersections.size(); i++) {
			Edge[] edges = is.intersectingEdges.get(i);

			Node c = (Node) is.intersections.get(i);
			QVector p = c.coord;

			QVector cp1 = (QVector) pointSet.get(edges[0].site1);
			QVector fp1 = edges[1].site1Node.coord;


			float[] d2 = new float[2];
			d2[0] = p.dist(cp1);
			d2[1] = p.dist(fp1);

			float dMin = findMin(d2);
			float dMax = findMax(d2);

			TwoCircle tc = new TwoCircle(p, dMin, dMax);
			if (tc.width <= minWidth || minWidth < 0) {
				minWidth = tc.width;
				best3 = tc;
			}
			if (tc.width >= maxWidth || maxWidth < 0) {
				maxWidth = tc.width;
				worst3 = tc;
			}
			cand.add(tc);
		}
		
		//System.out.println(best3);
		return cand;
	}

	private ArrayList detectCase2(FVoronoi fv, ArrayList pointSet) {
		if (fv == null || fv.fvNodes == null) {
			return null;
		}

		ArrayList cand = new ArrayList();
		double minWidth = -1, maxWidth = -1;

		for (int i = 0; i < fv.fvNodes.size(); i++) {
			Node n = (Node) fv.fvNodes.get(i);

			QVector p = new QVector();
			p.x = n.coord.x;
			p.y = n.coord.y;

			Iterator<Integer> iter = n.sites.iterator();
			int siteId = iter.next();
			QVector r = (QVector) fv.inputCH.get(siteId);

			float distMax = p.dist(r);
			float distMin = -1;

			for (int k = 0; k < pointSet.size(); k++) {
				QVector f = (QVector) pointSet.get(k);
				float dist = p.dist(f);
				if (dist < distMin || distMin < 0) {
					distMin = dist;
				}
			}

			TwoCircle tc = new TwoCircle(p, distMin, distMax);
			if (tc.width <= minWidth || minWidth < 0) {
				minWidth = tc.width;
				best2 = tc;
			}
			if (tc.width >= maxWidth || maxWidth < 0) {
				maxWidth = tc.width;
				worst2 = tc;
			}
			cand.add(tc);
		}

		return cand;
	}

	public ArrayList detectCase1(CVoronoi cv, ArrayList pointSet) {
		if (cv == null || cv.cvNodes == null) {
			return null;
		}

		ArrayList cand = new ArrayList();
		double minWidth = -1, maxWidth = -1;

		for (int i = 0; i < cv.cvNodes.size(); i++) {
			Node n = (Node) cv.cvNodes.get(i);

			QVector p = new QVector();
			p.x = n.coord.x;
			p.y = n.coord.y;

			Iterator<Integer> iter = n.sites.iterator();
			int siteId = iter.next();
			QVector r = (QVector) pointSet.get(siteId);

			float distMin = p.dist(r);
			float distMax = 0;

			for (int k = 0; k < pointSet.size(); k++) {
				QVector f = (QVector) pointSet.get(k);
				float dist = p.dist(f);
				if (dist > distMax) {
					distMax = dist;
				}
			}

			TwoCircle tc = new TwoCircle(p, distMin, distMax);
			if (tc.width <= minWidth || minWidth < 0) {
				minWidth = tc.width;
				best1 = tc;
			}
			if (tc.width >= maxWidth || maxWidth < 0) {
				maxWidth = tc.width;
				worst1 = tc;
			}
			cand.add(tc);
		}

		return cand;
	}
}