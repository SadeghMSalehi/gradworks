#!/usr/bin/python

from PySide.QtCore import *
from PySide.QtGui import *
from PySide import QtOpenGL

from visImageWindow import *

class VisImageWindow(QMainWindow, Ui_MainWindow):
	def __init__(self, parent=None):
		super(VisImageWindow, self).__init__(parent)
		self.setupUi(self)

if __name__ == "__main__":
	import sys
	_app = QApplication(sys.argv)
	_win = VisImageWindow()
	_win.show()
	sys.exit(_app.exec_())
