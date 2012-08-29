package com.intel.vpg;

public class Rect {
	public int x;
	public int y;
	public int w;
	public int h;

	public Rect(int[] data, int offset) {
		x = data[offset];
		y = data[offset + 1];
		w = data[offset + 2];
		h = data[offset + 3];
	}
	
	public Rect(Rect r) {
		x = r.x;
		y = r.y;
		w = r.w;
		h = r.h;
	}

	public void merge(Rect b) {
		x = Math.min(x, b.x);
		y = Math.min(y, b.y);
		w = Math.max(x + w, b.x + b.w) - x;
		h = Math.max(y + h, b.y + b.h) - y;
	}

	public void intersect(Rect b) {
		x = Math.max(x, b.x);
		y = Math.max(y, b.y);
		w = Math.min(x + w, b.x + b.w) - x;
		h = Math.min(y + h, b.y + b.h) - y;		
	}

	boolean valueInRange(int value, int min, int max) {
		return (value >= min) && (value <= max);
	}

	boolean overlap(Rect B) {
		boolean xOverlap = valueInRange(this.x, B.x, B.x + B.w)
				|| valueInRange(B.x, this.x, this.x + this.w);

		boolean yOverlap = valueInRange(this.y, B.y, B.y + B.h)
				|| valueInRange(B.y, this.y, this.y + B.h);

		return xOverlap && yOverlap;
	}

	public String toString() {
		return String.format("(%d, %d, %d, %d)", x, y, w, h);
	}

}
