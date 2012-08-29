package com.intel.vpg;

public class DoubleMatrix {
	public double[] data;
	public int nCols;
	public int nRows;

	public DoubleMatrix(int rows, int cols) {
		allocate(rows, cols);
	}

	public void allocate(int rows, int cols) {
		nRows = rows;
		nCols = cols;
		data = new double[rows * cols];
	}
}
