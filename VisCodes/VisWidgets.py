#!/usr/bin/env python

"""PySide port of the opengl/samplebuffers example from Qt v4.x"""

import sys
import math
from PySide.QtOpenGL import *
from PySide.QtCore import *
from PySide.QtGui import *

try:
	from OpenGL.GL import *
except ImportError:
	app = QApplication(sys.argv)
	QMessageBox.critical(None, "Opensamplebuffers",
													"PyOpenmust be installed to run this example.",
													QMessageBox.Ok | QMessageBox.Default,
													QMessageBox.NoButton)
	sys.exit(1)


class QProcessingBase(QWidget):
	GL_MULTISAMPLE = 0x809D

	def __init__(self, parent=None):
		#QGLWidget.__init__(self, QGLFormat(QGL.SampleBuffers), parent)
		QWidget.__init__(self, parent)
		self._smooth = 0
		self._timerId = 0
		self.grabKeyboard()
		self.grabMouse()
		self.setup()

	def initializeGL(self):
		pass

	def paintEvent(self,ev):
		p = QPainter(self)
		p.setRenderHint(QPainter.Antialiasing, self._smooth)
		self.draw(p)

	def keyReleaseEvent(sf,ev):
		sf.keyPressed(ev.text().encode('latin1'))

	def timerEvent(self, event):
		self.update()

class QProcessing(QProcessingBase):
	def frameRate(self, f):
		self._frameRate = f
		self.killTimer(self._timerId)
		self._timerId = self.startTimer(1000./self._frameRate)

	def smooth(self):
		self._smooth = 1
	
	def noSmooth(self):
		self._smooth = 0
	

class QImageViewer(QProcessing):
	def setup(self):
		self.frameNo = 0
		self.frameRate(30)
		self.smooth()

	def keyPressed(s,k):
		if k == 'q':
			sys.exit(0)

	def draw(s,p):
		s.frameNo += 1
		p.translate(s.width()/2.,s.height()/2.)
		for diameter in range(0,256,9):
			delta = abs((s.frameNo % 128) - diameter/2.)
			alpha = 255 - (delta*delta/4.) - diameter
			if alpha > 0:
				p.setPen(QPen(QColor(0,diameter/2.,127,alpha),3))
				p.drawEllipse(QRect(-diameter/2.,-diameter/2.,diameter,diameter))
