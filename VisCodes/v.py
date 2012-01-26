from PySide.QtCore import *
from PySide.QtGui import *
from PySide.QtOpenGL import *
from OpenGL import GL
import sys

class QProcessingBase(QGLWidget):
	def __init__(self, parent=None):
		QGLWidget.__init__(self,QGLFormat(QGL.SampleBuffers), parent)
		self._smooth = 0
		self.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding)
		self.startTimer(100)
		self.grabKeyboard()
		self.grabMouse()
		self.setup()

	def initializeGL(self):
		GL.glEnable(GL.GL_MULTISAMPLE)

	def animate(self):
		self.update()

	def keyReleaseEvent(sf,ev):
		sf.keyPressed(ev.text().encode('latin1'))

	def timerEvent(sf,ev):
		sf.update()

	def resizeGL(self,width,height):
		GL.glViewport(0,0,width,height)

	def paintEvent(self,ev):
		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glPushMatrix()
		colorPurple = QColor.fromCmykF(0.4,0.0,1.0,0.0)
		GL.glClearColor(.4,.0,1.,.0)
		GL.glShadeModel(GL.GL_SMOOTH)
		map(lambda x: GL.glEnable(x), \
			[GL.GL_DEPTH_TEST,GL.GL_LIGHTING,GL.GL_CULL_FACE,GL.GL_LIGHTING,GL.GL_LIGHT0,GL.GL_MULTISAMPLE])
		GL.glLightfv(GL.GL_LIGHT0, GL.GL_POSITION, (.5,5.,7.,1.))
		GL.glViewport(0,0,int(self.width()),int(self.height()))
		GL.glClear(GL.GL_COLOR_BUFFER_BIT|GL.GL_DEPTH_BUFFER_BIT)
		GL.glLoadIdentity()
		p = QPainter(self)
		p.begin(self)
		p.setRenderHint(QPainter.Antialiasing, self._smooth)

		s = self
		s.frameNo += 1
		p.translate(s.width()/2.,s.height()/2.)
		for diameter in range(0,256,9):
			delta = abs((s.frameNo % 128) - diameter/2.)
			alpha = 255 - (delta*delta/4.) - diameter
			if alpha > 0:
				p.setPen(QPen(QColor(0,diameter/2.,127,alpha),3))
				p.drawEllipse(QRect(-diameter/2.,-diameter/2.,diameter,diameter))
		p.end()

class QProcessing(QProcessingBase):
	def frameRate(self, f):
		self._frameRate = f
		self.startTimer(1000/self._frameRate)

	def smooth(self):
		self._smooth = 1
	
	def noSmooth(self):
		self._smooth = 0
	

class QImageViewer(QProcessing):
	def setup(self):
		self.frameRate(30)
		self.frameNo = 0

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
