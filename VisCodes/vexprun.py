#!/tools/ParaView-Build/bin/pvpython

# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'vexp.ui'
#
# Created: Wed Jan  4 22:31:20 2012
#      by: pyside-uic 0.2.11 running on PySide 1.0.9
#
# WARNING! All changes made in this file will be lost!

from PySide.QtGui import *
from PySide.QtCore import *

import matplotlib
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = "PySide"
import pylab

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import vexp
import vtk
import streamlines as sl

import numpy as np

class VEMainWindow(QMainWindow, vexp.Ui_MainWindow):
	def __init__(self, parent=None):
		super(VEMainWindow, self).__init__(parent)
		self.setupUi(self)

		self.ren = vtk.vtkRenderer()
		self.ren.SetBackground(0.2,0.2,0.2)

		self.i = self.qvtkWidget

		self.renWin = self.i.GetRenderWindow()
		self.renWin.SetAAFrames(5)
		self.renWin.AddRenderer(self.ren)
		self.renWin.Render()

		self.currLineId = 0

		self.polyGeom = sl.readVTP("/data/viscontest/Blade_Geometry.90.vtp")
		self.polyLine = None

		# generate the plot
		fig = Figure(figsize=(600,600), dpi=72, facecolor=(1,1,1), edgecolor=(0,0,0))
		self.ax = fig.add_subplot(111)
		self.ax.plot([0,1])
		# generate the canvas to display the plot
		canvas = FigureCanvas(fig)

		self.canvasWidget = canvas
		self.canvasWidget.setMinimumSize(QSize(400, 0))
		self.canvasWidget.setObjectName("canvasWidget")
		self.horizontalLayout_5.addWidget(self.canvasWidget)

	@Slot()
	def on_prevButton_clicked(self):
		if (self.currLineId > 0):
			self.currLineId = self.currLineId - 1
			self.ChangeLine()
			self.RenderLine()
	
	@Slot()
	def on_nextButton_clicked(self):
		if (self.currLineId < self.poly.GetNumberOfLines() - 1):
			self.currLineId = self.currLineId + 1
			self.ChangeLine()
			self.RenderLine()

	@Slot()
	def on_fitButton_clicked(self):
		"Fitting a sphere to the line hoping the sphere to capture a circular movement"
		ids = []

		lps = sl.firstLineAsArray(self.polyLine, ids)
		slps = sl.smoothLine(lps)
		dslps = (len(ids) - slps.shape[0]) / 2

		pts = self.polyLine.GetPoints()

		for (idx,j) in enumerate(range(dslps, len(lps)-dslps)):
			pts.SetPoint(ids[j], slps[idx,:])

	

	@Slot()
	def on_action_Save_triggered(self):
		if (self.polyLine is None):
			return

		(fileName, x) = QFileDialog.getSaveFileName(self, unicode("Save VTK File"), "/data/viscontest/DES_streamlines", unicode("VTK files (*.vtk)"))
		sl.writeVTK(fileName, self.polyLine)

	@Slot()
	def on_smoothButton_clicked(self):
		ids = []

		lps = sl.firstLineAsArray(self.polyLine, ids)
		slps = sl.smoothLine(lps)
		dslps = (len(ids) - slps.shape[0]) / 2

		pts = self.polyLine.GetPoints()

		for (idx,j) in enumerate(range(dslps, len(lps)-dslps)):
			pts.SetPoint(ids[j], slps[idx,:])

		self.RenderLine()

		return

	@Slot()
	def on_bPlotK_clicked(self):
		self.RenderLine()

	@Slot()
	def on_bPlotdK_clicked(self):
		self.RenderLine()

	@Slot()
	def on_bPlotKFFT_clicked(self):
		self.RenderLine()

	@Slot()
	def on_action_Open_triggered(self):
		(fileName, x) = QFileDialog.getOpenFileName(self, unicode("Open VTK File"), "/data/viscontest/DES_streamlines", unicode("VTK files (*.vtp *.vtk)"))
		self.poly = sl.readVTP(fileName)
		self.ibox.setHtml("# of lines: %d<br>" % (self.poly.GetNumberOfLines()))
		self.ChangeLine()
		self.RenderLine()

	def enumLine(self, poly, lineNo):
		iterLine = poly.GetLines()
		iterLine.InitTraversal()
		idList = vtk.vtkIdList()

		keepIter = 1
		lineCount = 0

		while (keepIter > 0):
			keepIter = iterLine.GetNextCell(idList)

			if (lineCount == lineNo):
				newLines = vtk.vtkCellArray()
				newLines.InsertNextCell(idList)
				return newLines

			lineCount = lineCount + 1
		return None

	def AddGeometry(self):
		m = vtk.vtkPolyDataMapper()
		m.SetInput(self.polyGeom)

		a = vtk.vtkActor()
		a.SetMapper(m)
		a.GetProperty().SetColor(0.8,0.8,0.8)
		a.GetProperty().SetOpacity(0.5)

		self.ren.AddActor(a)

	def ChangeLine(self):
		line = self.enumLine(self.poly, self.currLineId)

		p = vtk.vtkPolyData()
		p.SetPoints(self.poly.GetPoints())
		p.SetLines(line)

		self.polyLine = p


	def RenderLine(self):
		self.ren.RemoveAllViewProps()

		p = self.polyLine

		m = vtk.vtkPolyDataMapper()
		m.SetInput(p)

		a = vtk.vtkActor()
		a.SetMapper(m)

		self.AddGeometry()
		self.ren.AddActor(a)
		self.ren.ResetCamera()

		self.i.Render()

		ids = []
		line = sl.firstLineAsArray(p, ids)
		curvatures = sl.computeCurvatureForLine(line)


		self.ax.clear()

		if (self.bPlotK.isChecked()):
			self.ax.plot(curvatures)

		(p1, p2, dk) = sl.bracketPeaks(curvatures, -400, 400)

		if (self.bPlotdK.isChecked()):
			self.ax.plot(dk)

		if (self.bPlotKFFT.isChecked()):
			fftK = np.abs(np.fft.rfft(curvatures))
			self.ax.plot(fftK)

		self.canvasWidget.draw()


if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    _mainWindow = VEMainWindow()
    _mainWindow.show()
    sys.exit(app.exec_())

