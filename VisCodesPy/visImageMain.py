#!/usr/bin/python

from PySide.QtCore import *
from PySide.QtGui import *
from PySide import QtOpenGL

from visImageWindow import *
from QProcessing import *

class CircleTest(QProcessing):
	def setup(self):
		self.frameNo = 0
		self.frameRate(30)
		self.smooth()

	def keyPressed(s,k):
		if k == 'q':
			sys.exit(0)

	def draw(s):
		s.background(255)
		s.frameNo += 1
		s.translate(s.width()/2.,s.height()/2.)
		"""
		for diameter in range(0,256,9):
			delta = abs((s.frameNo % 128) - diameter/2.)
			alpha = 255 - (delta*delta/4.) - diameter
			if alpha > 0:
				s.strokeWeight(3)
				s.stroke(0,diameter/2.,127,alpha)
				s.ellipse(-diameter/2.,-diameter/2.,diameter,diameter)
		s.stroke(0)
		s.line(30, 20, 85, 20);
		s.stroke(126);
		s.line(85, 20, 85, 75);
		s.stroke(255);
		s.line(85, 75, 30, 75);
		"""
		s.ellipseMode(CENTER)
		s.arc(50, 55, 50, 50, 0, PI/2);
		s.noFill();
		s.arc(50, 55, 60, 60, PI/2, PI);
		s.arc(50, 55, 70, 70, PI, TWO_PI-PI/2);
		s.arc(50, 55, 80, 80, TWO_PI-PI/2, TWO_PI);


class VisImageWindow(QMainWindow, Ui_MainWindow):
	def __init__(self, parent=None):
		super(VisImageWindow, self).__init__(parent)
		self.setupUi(self)
		self.widget.setProcessing(CircleTest(self.widget))
	
	@Slot()
	def on_action_Close_triggered(self):
		sys.exit(1)

if __name__ == "__main__":
	import sys
	_app = QApplication(sys.argv)
	_win = VisImageWindow()
	_win.show()
	sys.exit(_app.exec_())
