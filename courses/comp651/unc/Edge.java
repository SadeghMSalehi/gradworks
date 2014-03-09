package unc;

public class Edge {
	public Node origin;
	public Node dest;
	public int site1, site2;
	public Node site1Node;
	public Node site2Node;

	// public QVector dir;

	public Edge(Node o, Node d) {
		origin = o;
		dest = d;
	}

	public Edge(Node o, Node d, int site1, int site2) {
		origin = o;
		dest = d;
		this.site1 = site1;
		this.site2 = site2;
	}
	
	public Edge(float x1, float y1, float x2, float y2) {
		origin = new Node(x1, y1);
		dest = new Node(x2, y2);
	}

	public Edge(QVector p, QVector q) {
		origin = new Node(p);
		dest = new Node(q);
	}

	public Edge(double x1, double y1, double x2, double y2) {
		this((float) x1, (float) y1, (float) x2, (float) y2);
	}

	//
	// public void dir(QVector o, QVector d) {
	// this.dir = QVector.sub(d, o);
	// this.dir.normalize();
	// }
	//
	// public QVector intersect(Edge o) {
	// float x1 = origin.coord.x, y1 = origin.coord.y, x2 = dest.coord.x, y2 =
	// dest.coord.y;
	// float x3 = o.origin.coord.x, y3 = o.origin.coord.y, x4 = o.dest.coord.x,
	// y4 = o.dest.coord.y;
	//
	// double ccw = (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4);
	// if (ccw == 0) {
	// System.out.println("Parallel:" + ccw);
	// return null;
	// }
	//
	// double px = ((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/ccw;
	// double py = ((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/ccw;
	//
	// boolean isOn1 = Math.min(x1, x2)-TOL <= px && px <= Math.max(x1, x2)+TOL
	// && Math.min(y1, y2)-TOL <= py && py <= Math.max(y1, y2)+TOL;
	// boolean isOn2 = Math.min(x3, x4)-TOL <= px && px <= Math.max(x3, x4)+TOL
	// && Math.min(y3, y4)-TOL <= py && py <= Math.max(y3, y4)+TOL;
	// if (!isOn1 || !isOn2) {
	// System.out.println(px + "/" + py);
	// return null;
	// }
	//
	// QVector pOut = new QVector();
	// pOut.x = (float) px;
	// pOut.y = (float) py;
	//
	// return pOut;
	// }
	//
	// public boolean isOn(float x, float y) {
	// float x1 = origin.coord.x, y1 = origin.coord.y, x2 = dest.coord.x, y2 =
	// dest.coord.y;
	// double ccw = (x1-x2)*(y2-y)-(y1-y2)*(x2-x);
	// return Math.abs(ccw) < TOL && Math.min(x1, x2) <= x && x <= Math.max(x1,
	// x2) && Math.min(y1, y2) <= y && y <= Math.max(y1, y2);
	// }
	//
	// /*
	// * the edge is turning left (+), turning right (-), or parallel (0)
	// */
	// public double ccw(Edge o) {
	// double x1 = origin.coord.x, y1 = origin.coord.y, x2 = dest.coord.x, y2 =
	// dest.coord.y;
	// double x3 = o.origin.coord.x, y3 = o.origin.coord.y, x4 =
	// o.origin.coord.x, y4 = o.origin.coord.y;
	//
	// double ccw = (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4);
	// return ccw;
	// }
	//
	public Edge next;
	public Edge prev;
	// public Edge sym;
	// public boolean infOut;
	// public boolean infIn;
}