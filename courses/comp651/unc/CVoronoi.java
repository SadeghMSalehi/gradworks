package unc;
import java.util.ArrayList;
import java.util.List;

import be.humphreys.simplevoronoi.GraphEdge;
import be.humphreys.simplevoronoi.Voronoi;

public class CVoronoi {
	ArrayList pointSet;
	ArrayList cvEdges;
	ArrayList cvNodes;
	List<double[]> nodes = new ArrayList<double[]>();
	float xmin, ymin, xmax, ymax;

	public CVoronoi(ArrayList pointSet, float x0, float y0, float x1, float y1) {
		this.pointSet = pointSet;
		xmin = x0;
		ymin = y0;
		xmax = x1;
		ymax = y1;

		cvEdges = new ArrayList();
		cvNodes = new ArrayList();

		Voronoi v = new Voronoi(0.1);
		double[] xValuesIn = new double[pointSet.size()];
		double[] yValuesIn = new double[pointSet.size()];

		for (int i = 0; i < pointSet.size(); i++) {
			QVector s = (QVector) pointSet.get(i);
			xValuesIn[i] = s.x;
			yValuesIn[i] = s.y;
		}

		if (xValuesIn.length > 0) {
			List<GraphEdge> l = v.generateVoronoi(xValuesIn, yValuesIn, x0, x1,
					y0, y1);

			for (GraphEdge e : l) {
				Edge m = new Edge(e.x1, e.y1, e.x2, e.y2);
				m.site1 = e.site1;
				m.site2 = e.site2;
				cvEdges.add(m);

				double px = 0, py = 0;
				int foundNodeIndex = -1;
				if (e.x1 == xmin || e.x1 == xmax || e.y1 == ymin || e.y1 == ymax) {
					if (e.x2 == xmin || e.x2 == xmax || e.y2 == ymin || e.y2 == ymax) {
						continue;
					} else {
						foundNodeIndex = findNode(e.x2, e.y2);
						px = e.x2;
						py = e.y2;
					}
				} else {
					foundNodeIndex = findNode(e.x1, e.y1);
					px = e.x1;
					py = e.y1;
				}
				
				if (foundNodeIndex < 0) {
					nodes.add(new double[] { px, py });
					Node n = new Node((float) px, (float) py);
					n.sites.add(e.site1);
					n.sites.add(e.site2);
					cvNodes.add(n);
				} else {
					Node n = (Node) cvNodes.get(foundNodeIndex);
					n.sites.add(e.site1);
					n.sites.add(e.site2);
				}

			}
		}
	}

	private boolean compare(double x, double y) {
		return (Math.abs(x-y) < 1e-6);
	}
	private int findNode(double x1, double y1) {
		int foundNodeIndex = 0;
		for (double[] n : nodes) {
			boolean found = compare(n[0], x1) && compare(n[1],y1);
//			System.out.println(x1 + "," + y1 + " : " + n[0] + "," + n[1] + " = " + (found ? "match" : "mismatch"));
			if (found) {
				return foundNodeIndex;
			}
			foundNodeIndex++;
		}
		return -1;
	}
}