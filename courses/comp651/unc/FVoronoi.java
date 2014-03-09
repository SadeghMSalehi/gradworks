package unc;

import java.util.ArrayList;
import java.util.Comparator;

import com.vividsolutions.jts.algorithm.LineIntersector;
import com.vividsolutions.jts.algorithm.NonRobustLineIntersector;
import com.vividsolutions.jts.geom.Coordinate;

public class FVoronoi {
	public static float TOL = 1e-6f;

	public ArrayList inputCH;
	public ArrayList pointSet;

	public int nSites;
	public Node[] sites;
	public Node[] anchors;
	public float x0, y0, x1, y1;
	public ArrayList fvNodes, fvEdges;

	public QVector infinityAnchor(QVector p, QVector r) {
		QVector pr = QVector.sub(r, p);
		QVector npr = new QVector();

		npr.x = -pr.y * ((float) Math.max(x1 - x0, y1 - y0));
		npr.y = pr.x * ((float) Math.max(x1 - x0, y1 - y0));

		QVector mpr = new QVector();
		mpr.x = (p.x + r.x) / 2;
		mpr.y = (p.y + r.y) / 2;

		// npr.scaleTo((float) 2*Math.max(x1-x0, y1-y0));
		QVector tnpr = QVector.add(mpr, npr);

		Edge ep = new Edge(mpr, tnpr);

		Edge e1 = new Edge(x0, y0, x1, y0);
		Edge e2 = new Edge(x1, y0, x1, y1);
		Edge e3 = new Edge(x1, y1, x0, y1);
		Edge e4 = new Edge(x0, y1, x0, y0);

		QVector vp1 = intersectEdge(e1, ep);
		QVector vp2 = intersectEdge(e2, ep);
		QVector vp3 = intersectEdge(e3, ep);
		QVector vp4 = intersectEdge(e4, ep);

		QVector pClip = null;
		if (vp1 != null) {
			pClip = vp1;
		} else if (vp2 != null) {
			pClip = vp2;
		} else if (vp3 != null) {
			pClip = vp3;
		} else if (vp4 != null) {
			pClip = vp4;
		} else {
			System.out.println("Non-existing VP");
		}

		return pClip;
	}

	private QVector intersectEdge(Edge e1, Edge ep) {
		Coordinate p1 = new Coordinate();
		Coordinate p2 = new Coordinate();
		Coordinate p3 = new Coordinate();
		Coordinate p4 = new Coordinate();

		p1.x = e1.origin.coord.x;
		p1.y = e1.origin.coord.y;
		p2.x = e1.dest.coord.x;
		p2.y = e1.dest.coord.y;
		p3.x = ep.origin.coord.x;
		p3.y = ep.origin.coord.y;
		p4.x = ep.dest.coord.x;
		p4.y = ep.dest.coord.y;

		LineIntersector li = new NonRobustLineIntersector();
		li.computeIntersection(p1, p2, p3, p4);

		if (li.hasIntersection()) {
			Coordinate p = li.getIntersection(0);
			QVector pOut = new QVector();
			pOut.x = (float) p.x;
			pOut.y = (float) p.y;
			return pOut;
		}

		return null;
	}

	public FVoronoi(ArrayList pointSet, float x0, float y0, float x1, float y1) {
		if (pointSet.size() < 3) {
			return;
		}

		this.pointSet = pointSet;
		this.inputCH = constructConvexHull();
		// this.inputCH.remove(chull.size()-1);

		fvNodes = new ArrayList();
		fvEdges = new ArrayList();

		nSites = inputCH.size() - 1;
		if (nSites < 3) {
			return;
		}

		this.x0 = x0;
		this.y0 = y0;
		this.x1 = x1;
		this.y1 = y1;

		sites = new Node[nSites];
		anchors = new Node[nSites];

		for (int i = 0; i < nSites; i++) {
			sites[i] = new Node((QVector) inputCH.get(i));
			sites[i].id = i;

			if (i > 0) {
				sites[i].prev = sites[i - 1];
				sites[i - 1].next = sites[i];
			}
		}
		sites[0].prev = sites[nSites - 1];
		sites[nSites - 1].next = sites[0];

		initAnchors();

		int n = nSites;
		Node head = sites[0];

		do {
			// find maximum circle
			Node p = null;
			QVector maxC = null;

			Node iter = head;

			double maxR = 0;
			double maxA = 0;

			// find p maximizing radius and angle
			do {
				QVector p1 = iter.coord;
				QVector p2 = iter.next.coord;
				QVector p3 = iter.next.next.coord;

				QVector c = computeCircle(p1, p2, p3);
				double angle = QVector.angleBetween(QVector.sub(p1, p2),
						QVector.sub(p3, p2));

				if (c == null) {
					// three co-circular point
					continue;
				}
				if (c.z > maxR) {
					maxR = c.z;
					maxA = angle;

					maxC = c;
					maxC.z = 0;
					p = iter.next;
				} else if (c.z == maxR) {
					if (maxA < angle) {
						maxA = angle;
						maxC = c;
						maxC.z = 0;
						p = iter.next;
					}
				}
				iter = iter.next;
			} while (iter != head);

			Node q = p.prev;
			Node c = new Node(maxC);
			c.prev = q;
			c.next = q.next;

			Node vp = anchors[p.id];
			Node vq = anchors[q.id];

			fvNodes.add(c);
			c.sites.add(q.id);
			c.sites.add(p.id);
			c.sites.add(p.next.id);

			Edge e1 = new Edge(c, vp);
			e1.site1Node = p;
			e1.site2Node = p.next;
			fvEdges.add(e1);

			Edge e2 = new Edge(c, vq);
			e2.site1Node = q;
			e2.site2Node = q.next;
			fvEdges.add(e2);

			anchors[q.id] = c;

			// remove p
			q.next = p.next;
			q.next.prev = q;

			if (p == head) {
				head = q;
			}

			n = n - 1;
		} while (n > 2);

		Edge e = new Edge(anchors[head.id], anchors[head.next.id]);
		e.site1Node = head;
		e.site2Node = head.next;
		fvEdges.add(e);
	}

	// be careful when to use the returning QVector. It has z-component
	public QVector computeCircle(QVector p1, QVector p2, QVector p3) {
		double det = (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y)
				* (p3.x - p1.x);
		if (det != 0) {
			double rx = (p2.x * p2.x - p1.x * p1.x + p2.y * p2.y - p1.y * p1.y) / 2;
			double ry = (p3.x * p3.x - p1.x * p1.x + p3.y * p3.y - p1.y * p1.y) / 2;

			QVector c = new QVector();
			c.x = (float) (((p3.y - p1.y) * rx - (p2.y - p1.y) * ry) / det);
			c.y = (float) ((-(p3.x - p1.x) * rx + (p2.x - p1.x) * ry) / det);
			c.z = 0;

			assert (p1.z == 0);
			c.z = c.dist(p1);

			return c;
		}
		return null;
	}

	private void initAnchors() {
		for (int i = 0; i < sites.length; i++) {
			Node p = sites[i];
			Node r = sites[(i + 1) % nSites];

			QVector pr = infinityAnchor(p.coord, r.coord);
			anchors[i] = new Node(pr);
			anchors[i].prev = p;
			anchors[i].next = r;
		}

	}

	public double ccw(QVector p1, QVector p2, QVector p3) {
		return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
	}

	/**
	 * find the most left bottom point find the second left-most point and
	 * connect, compare the rest of
	 */
	public ArrayList constructConvexHull() {
		ArrayList convexPoints = new ArrayList();

		if (pointSet.size() < 2) {
			return null;
		}

		ArrayList convexHull = (ArrayList) pointSet.clone();

		Comparator zcomp = new Comparator() {
			public int compare(Object o1, Object o2) {
				QVector p1 = (QVector) o1;
				QVector p2 = (QVector) o2;

				if (p1.z < p2.z) {
					return -1;
				} else if (p1.z == p2.z) {
					return 0;
				} else {
					return 1;
				}
			}
		};

		// find left-bottom node
		int leftBottom = 0;
		for (int i = 0; i < convexHull.size(); i++) {
			QVector p1 = (QVector) convexHull.get(leftBottom);
			QVector p2 = (QVector) convexHull.get(i);

			if (p1.y > p2.y) {
				leftBottom = i;
			} else if (p1.y == p2.y) {
				if (p1.x > p2.x) {
					leftBottom = i;
				}
			}
		}

		// compute angles
		for (int i = 0; i < convexHull.size(); i++) {
			if (i == leftBottom) {
				continue;
			}
			QVector p0 = (QVector) convexHull.get(leftBottom);
			QVector pi = (QVector) convexHull.get(i);

			double x2x1 = pi.x - p0.x;
			double y2y1 = pi.y - p0.y;
			double dist = Math.sqrt(x2x1 * x2x1 + y2y1 * y2y1);
			double cos = x2x1 / dist;
			pi.z = (float) (1 - cos) + 1;
		}

		convexHull = sortPoints(convexHull, zcomp);
		convexHull.add(convexHull.get(0));

		// graham's scan
		int m = 1;
		for (int i = 2; i < convexHull.size(); i++) {
			QVector p0 = (QVector) convexHull.get(m - 1);
			QVector p1 = (QVector) convexHull.get(m);
			QVector pi = (QVector) convexHull.get(i);

			while (ccw(p0, p1, pi) < 0) {
				// collinear with bottom-left node
				if (m == 1) {
					convexHull.set(m, pi);
					convexHull.set(i, p1);

					i += 1;
				} else {
					m -= 1;
				}

				p0 = (QVector) convexHull.get(m - 1);
				p1 = (QVector) convexHull.get(m);
				pi = (QVector) convexHull.get(i);
			}

			m += 1;
			p1 = (QVector) convexHull.get(m);
			convexHull.set(m, pi);
			convexHull.set(i, p1);
		}

		for (int i = 0; i <= m; i++) {
			convexPoints.add(convexHull.get(i));
		}

		// System.out.println(convexPoints);
		return convexPoints;
	}

	public ArrayList sortPoints(ArrayList ps, Comparator comp) {
		int len = ps.size();
		if (len < 2) {
			return ps;
		}
		if (len <= 2) {
			QVector p1 = (QVector) ps.get(0);
			QVector p2 = (QVector) ps.get(1);

			if (comp.compare(p1, p2) < 0) {
				return ps;
			} else {
				ps.set(0, p2);
				ps.set(1, p1);
				return ps;
			}
		}

		int pivot = (int) (Math.random() * len);

		QVector pivotPoint = (QVector) ps.get(pivot);
		ArrayList small = new ArrayList();
		ArrayList large = new ArrayList();

		for (int i = 0; i < ps.size(); i++) {
			if (i == pivot) {
				continue;
			}
			QVector p = (QVector) ps.get(i);
			if (comp.compare(p, pivotPoint) < 0) {
				small.add(p);
			} else {
				large.add(p);
			}
		}

		ArrayList sortedSmall = sortPoints(small, comp);
		ArrayList sortedLarge = sortPoints(large, comp);

		sortedSmall.add(pivotPoint);
		sortedSmall.addAll(sortedLarge);

		return sortedSmall;
	}

	public ArrayList shufflePoints(ArrayList ps) {
		for (int i = 0; i < ps.size(); i++) {
			int r = (int) (Math.random() * (ps.size() - i));
			QVector pi = (QVector) ps.get(i);
			QVector pr = (QVector) ps.get(r);
			ps.set(i, pr);
			ps.set(r, pi);
		}
		return ps;
	}
}
