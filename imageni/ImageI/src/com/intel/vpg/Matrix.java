package com.intel.vpg;

public class Matrix {
	public float[] data;
	private int nCols;
	private int nRows;

	public Matrix(int rows, int cols) {
		allocate(rows, cols);
	}

	public void allocate(int rows, int cols) {
		setNumRows(rows);
		setNumCols(cols);
		data = new float[rows * cols];
	}

	public int getNumRows() {
		return nRows;
	}

	public void setNumRows(int nRows) {
		this.nRows = nRows;
	}

	public int getNumCols() {
		return nCols;
	}

	public void setNumCols(int nCols) {
		this.nCols = nCols;
	}
}
