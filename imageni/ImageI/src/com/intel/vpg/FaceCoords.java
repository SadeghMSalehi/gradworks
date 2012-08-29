package com.intel.vpg;

import java.util.ArrayList;

public class FaceCoords {
	private Rect[] _data;

	public FaceCoords(int[] data) {
		init(data);
	}

	private void init(int[] data) {
		ArrayList<Rect> tmp = new ArrayList<Rect>();
		for (int i = 0; i < data.length; i += 4) {
			Rect r = new Rect(data, i);
			if (tmp.size() == 0) {
				tmp.add(r);
			} else {
				boolean overlap = false;
				for (int j = 0; j < tmp.size(); j++) {
					if (tmp.get(j).overlap(r)) {
						tmp.get(j).merge(r);
						overlap = true;
					}
				}
				if (!overlap) {
					tmp.add(r);
				}
			}
		}
		_data = tmp.toArray(new Rect[tmp.size()]);
	}

	public int size() {
		return _data.length;
	}
	
	public Rect getRect(int i) {
		return _data[i];
	}
}
