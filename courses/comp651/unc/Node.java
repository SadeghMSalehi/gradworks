package unc;
import java.util.HashSet;

public class Node {
	public QVector coord;
	public Node prev;
	public Node next;

	public int id;
	public HashSet<Integer> sites;
	
	public Node() {
		sites = new HashSet<Integer>();
	};

	public Node(QVector coord) {
		this();
		this.coord = coord;
	}

	public Node(float x, float y) {
		this();
		this.coord = new QVector();
		this.coord.x = x;
		this.coord.y = y;
	}
}