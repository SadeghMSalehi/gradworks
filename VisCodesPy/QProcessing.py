#!/usr/bin/env python

"""PySide port of the opengl/samplebuffers example from Qt v4.x"""

import sys
from PySide.QtOpenGL import *
from PySide.QtCore import *
from PySide.QtGui import *
from math import *

try:
	from OpenGL.GL import *
except ImportError:
	print "OpenGL cannot be found"
	sys.exit(1)


"""
Proxy class for QProcessing
Accept QT events and hand over to user's Processing class
"""
class QProcessingWidget(QWidget):
	GL_MULTISAMPLE = 0x809D

	def __init__(self, parent=None):
		#QGLWidget.__init__(self, QGLFormat(QGL.SampleBuffers), parent)
		QWidget.__init__(self, parent)
		self._timerId = 0
		self.setFocusPolicy(Qt.StrongFocus)
		self.setFocus()
		self.P = None

	def setProcessing(self,p):
		self.P = p

	def initializeGL(self):
		pass
	
	def resizeEvent(self,ev):
		self.P.resized()

	def paintEvent(self,ev):
		self.P.draw()
		painter = QPainter(self)
		pixmap = self.P.hold()
		painter.drawPixmap((self.width()-pixmap.width())/2,0,pixmap)
		painter.end()
		self.P.resume()

	def keyReleaseEvent(self,ev):
		self.P.keyPressed(ev.text().encode('latin1'))

	def timerEvent(self, event):
		self.update()


PI = pi
TWO_PI = 2*pi
CORNER = 0
CENTER = 1
MITER = 2
BEVEL = 3
ROUND = 4


"""
Actual Processing class
"""
class QProcessing(QObject):
	def __init__(self,widget):
		self.W = widget
		self._smooth = 0
		self._bgColor = QColor.fromRgb(255,255,255)
		self._stroke = (0)
		self._strokeWeight = 1
		self._pixmap = QPixmap(512,512)
		self._translateX = 0
		self._translateY = 0
		self._timerId = 1
		self._ellipseMode = CORNER
		self.P = QPainter(self._pixmap)
		self.P.setRenderHint(QPainter.Antialiasing,self._smooth)
		self.P.setRenderHint(QPainter.TextAntialiasing,self._smooth)
		self.setup() 
		
	def hold(self):
		self.P.end()
		return self._pixmap

	def resume(self):
		self.P.begin(self._pixmap)
		self.P.setRenderHint(QPainter.Antialiasing,self._smooth)
		self.P.setRenderHint(QPainter.TextAntialiasing,self._smooth)

	def resized(self):
		self.P.end()
		self._pixmap = QPixmap(self.width(),self.height())
		self.P.begin(self._pixmap)

	def width(self):
		return self.W.width()

	def height(self):
		return self.W.height()

	def background(self,*args):
		if len(args) == 1:
			self._bgColor = QColor.fromRgb(args[0],args[0],args[0])
		elif len(args) == 3:
			self._bgColor = QColor.fromRgb(args[0],args[1],args[2])
		elif len(args) == 4:
			self._bgColor = QColor.fromRgb(args[0],args[1],args[2],args[3])
		self.P.fillRect(0,0,self.width(),self.height(),self._bgColor)
	
	def smooth(self):
		self._smooth = 1
		self.P.setRenderHint(QPainter.Antialiasing,self._smooth)
		self.P.setRenderHint(QPainter.TextAntialiasing,self._smooth)

	def noSmooth(self):
		self._smooth = 0
		self.P.setRenderHint(QPainter.Antialiasing,self._smooth)
		self.P.setRenderHint(QPainter.TextAntialiasing,self._smooth)

	def ellipseMode(self,mode):
		if mode == CENTER:
			self._ellipseMode = mode
		else:
			self._ellipseMode = CORNER
	
	def ellipse(self,x,y,w,h):
		if self._ellipseMode == CENTER:
			self.P.drawEllipse(QPoint(x,y),w,h)
		else:
			self.P.drawEllipse(x,y,w,h)

	def stroke(self,*args):
		p = self.P
		newPen = QPen(self.P.pen())
		if len(args) == 1:
			newPen.setColor(QColor.fromRgb(args[0]))
		elif len(args) == 3:
			newPen.setColor(QColor.fromRgb(args[0],args[1],args[2]))
		elif len(args) >= 4:
			newPen.setColor(QColor.fromRgb(args[0],args[1],args[2],args[3]))
		p.setPen(newPen)

	def strokeJoin(self,w):
		pass

	def strokeWeight(self,w):
		newPen = QPen()
		newPen.setBrush(QBrush(self.P.pen().color()))
		newPen.setWidth(w)
		self.P.setPen(newPen)

	def frameRate(self, f):
		self.W.killTimer(self._timerId)
		self._frameRate = f
		self._timerId = self.W.startTimer(1000./self._frameRate)

	def translate(self,x,y):
		self._translateX = x
		self._translateY = y
		self.P.translate(x,y)

 	def noLoop(self):
		self.W.killTimer(self._timerId)

	def line(self,*p):
		if (len(p) == 4):
			self.P.drawLine(p[0],p[1],p[2],p[3])
		elif len(p) == 6:
			pass
	
	def point(self,*p):
		if (len(p) == 2):
			self.P.drawPoint(p[0],p[1])

	def arc(self,*p):
		if (len(p) == 6):
			s = (p[4]-pi/2)*180*16./pi
			e = (p[5]-p[4])*180*16./pi
			w = p[2]
			h = p[3]
			if (self._ellipseMode == CENTER):
				x = p[0]-w/2
				y = p[1]-h/2
			else:
				x = p[0]
				y = p[1]
		self.P.drawArc(x,y,w,h,s,e)

	def noFill(self):
		pass
