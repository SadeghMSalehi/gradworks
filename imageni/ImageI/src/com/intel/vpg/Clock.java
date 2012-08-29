package com.intel.vpg;

public class Clock {
	private long _startTime;
	
	public Clock() {
		tick();
	}

	public void tick() {
		_startTime = System.currentTimeMillis();
	}

	public int tock() {
		return (int) (System.currentTimeMillis() - _startTime);
	}
}
