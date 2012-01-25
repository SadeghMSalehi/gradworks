package unc;
import java.util.ArrayList;

import com.vividsolutions.jts.algorithm.LineIntersector;
import com.vividsolutions.jts.algorithm.NonRobustLineIntersector;
import com.vividsolutions.jts.geom.Coordinate;


public class Intersector {

	public ArrayList intersections;
	public ArrayList<Edge[]> intersectingEdges = new ArrayList<Edge[]>();
	
	public Intersector(ArrayList cvEdges, ArrayList fvEdges) {
		intersections = new ArrayList();
		
		for (int i = 0; i < cvEdges.size(); i++) {
			Edge c = (Edge) cvEdges.get(i);
			for (int j = 0; j < fvEdges.size(); j++) {
				Edge f = (Edge) fvEdges.get(j);
				
				QVector p = intersectEdge(c, f);
				if (p != null) {
					Node n = new Node(p);
					intersections.add(n);
					intersectingEdges.add(new Edge[] { c, f } );
				}
			}
		}
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
	
}
