#!/usr/bin/python

from pyprocessing import *
import scipy.optimize as op
import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from optparse import OptionParser
import Tkinter as tk
import tkMessageBox

Xg = None
Xs = None
Xr = None
Xi = None
X0 = None
Xe = np.asarray((1,0))
Xc = None
rt = 0
Theta = 0
Axis = [0,0]

def setup():
	size(500,500)
	smooth()
	ellipseMode(CENTER)
	global Xg
	Xg = np.genfromtxt("point.txt")
	ellipseFit(Xg)

def rotation(x,y):
	return atan2(y[1],y[0])	

def mouseDragged():
	global Xe
	Xe[0] = (mouse.x - width/2)
	Xe[1] = (mouse.y - height/2)
	n = np.linalg.norm(Xe)
	Xe = Xe / n
	
def rotatingEllipse():
	translate(width/2,height/2)
	r = rotation(None,Xe)
	rotate(rotation(None,Xe))
	ellipse(0,0,50,30)
	translate(0,0)
	rotate(0)
	return

def draw():
	background(255)
	noStroke()
	fill(0xFF0000FF)
	if Xg is not None:
		if Xc is not None:
			fill(0xFF0000FF)
			noFill()
			stroke(0)
			strokeWeight(1)
			translate(Xc[0],Xc[1])
			rotate(Theta)
			ellipse(0,0,2*Axis[0],2*Axis[1])
			rotate(-Theta)
			translate(-Xc[0],-Xc[1])
			noStroke()
		#stroke(0,0,0)
		#strokeWeight(1)
		fill(0xFF7FB7C7)
		for p in Xg:
			ellipse(p[0],p[1],7,7)
		if Xs is not None:
			noFill()
			#stroke(0xFF962918)
			#strokeWeight(3)
			ellipse(Xs[0],Xs[1],2*Xs[2],2*Xs[2])
		if Xr is not None:
			#noFill()
			#stroke(0xFFA60000)
			#strokeWeight(5)
			fill(0x66A60000)
			for p in Xr:
				ellipse(p[0],p[1],9,9)
		fill(0)
		if Axis is None:
			text("#:%d" % (Xg.shape[0]), 10, 10 + textAscent())
		else:
			text("#:%d fit: radius(%f,%f)" % (Xg.shape[0],Axis[0],Axis[1]), 10, 10 + textAscent())
	
def f(X):
	global Xi
	return np.power(np.sqrt(((Xi-X[0:2])*(Xi-X[0:2])).sum(1)) - X[2], 2).sum(0)

def randomSample(data,K):
	"sample K number of rows from data"
	N = data.shape[0]
	return data[np.random.random_integers(0,N-1,K),:]
	
def asEquation(C):
	r = ""
	base = ["x^2","x*y","y^2","x","y","1"]
	for (c,b) in zip(C,base):
		if (r == ""):
			r = "%f*%s" % (c,b)
		elif b == "1":
			r = r + "%+f" % (c)
		elif np.abs(c) <= 1e-15:
			continue
		else:
			r = r + "%+f*%s" % (c,b)
	return r

def RotationOfAxes(C):
	(a,b,c,d,e,f) = (C[0],C[1],C[2],C[3],C[4],C[5])
	t = atan2(b,a-c)/2.0 if (a!=c) else math.pi/4.0
	ct = math.cos(t)
	st = math.sin(t)
	(ct2,st2,ctst) = (ct*ct,st*st,ct*st)
	eq = "%f*x^2+%f*x*y+%f*y^2+%f*x+%f*y+%f" % (a,b,c,d,e,f)
	coeff = (a*ct2+b*ctst+c*st2,b*(ct2-st2)-2*(a-c)*ctst,a*st2-b*ctst+c*ct2,d*ct+e*st,e*ct-d*st,f)
	#neq = "%f*x^2+%f*x*y+%f*y^2+%f*x+%f*y+%f" % (A,B,C,D,E,F)
	print asEquation(coeff)
	return coeff

def ComputeEllipseParameter(C):
	(a11,a12,a22,b1,b2,c) = (C[0],C[1]/2,C[2],C[3],C[4],C[5])
	"Center K"
	(k1,k2) = np.asarray((a22*b1-a12*b2,a11*b2-a12*b1))/(2*(a12*a12-a11*a22))
	m = (a11*k1*k1+2*a12*k1*k2+a22*k2*k2-c)
	"Quadric form M"
	(m11,m12,m22)=(a11/m,a12/m,a22/m)
	""" direct computation
	"Eigenvalues of M"
	d = np.sqrt((m11-m22)*(m11-m22)+4*m12*m12)
	(l1,l2) = (((m11+m22)+d)/2.0,((m11+m22-d)/2.0))
	"Semi-major and minor axis"
	(a,b) = 1/np.sqrt((l2,l1))
	"""
	w,U = np.linalg.eig(np.linalg.inv(((m11,m12),(m12,m22))))
	(a,b) = np.sqrt(w)
	r = math.atan2(U[1,0],U[0,0])
	return ((a,b),(k1,k2),r)

def ellipseCostFunc(C,(x,y)):
	return C[0]*x*x+C[1]*x*y+C[2]*y*y+C[3]*x+C[4]*y+C[5]

def ellipseFit(data):
	global Xc,Theta,Axis
	if (data.shape[0] > 10):
		D = np.zeros((data.shape[0],6),dtype=np.float64)
		for i,(x,y) in enumerate(data):
			D[i,:] = [x*x,x*y,y*y,x,y,1]
		S = np.dot(D.transpose(),D)
		Sinv = np.linalg.inv(S)
		C = np.mat([[0,0,2,0,0,0],[0,-1,0,0,0,0],[2,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])
		A = np.dot(Sinv,C)
		(w,U) = np.linalg.eig(A)
		r = np.nonzero(w > 0)[0][0]
		u = U[:,r]
		uCu = np.dot(u.transpose(),np.dot(C,u))[0,0]
		a = np.asarray(u.transpose()/sqrt(uCu))[0]
		(Axis,Xc,Theta) = ComputeEllipseParameter(a)
		"return whose constant term is -1"
		return -a/a[5]
	else:
		print "I need at least 10 number of data points"
		
def findInliers(data, sample, model):
		d = sample-model[0:2]
		d = np.abs(np.sqrt((d*d).sum(1))-model[2])
		r = np.median(d)
		print r
		D = data-model[0:2]
		D = np.abs(np.sqrt((D*D).sum(1))-model[2])
		print D
		return np.nonzero(D<r)
	
def keyPressed():
	global Xg,Xr,Axis,Xc,Theta
	if key.char in "sS":
		if Xg is not None:
			np.savetxt("point.txt", Xg)
	elif key.char in "cC":
		print "Screen captured at screen.png"
		save("screen.png")
	elif key.char in "lL":
		print "Loading previously saved point.txt"
		Xg = np.genfromtxt("point.txt")
	elif key.char in "oO":
		print "run optimization"
		optimize()
	elif key.char in "hH":
		print "(S)ave/(C)apture/(L)oad/(O)ptimize/(H)elp"
	elif key.char in "nN":
		print "New"
		Xg = None
		Xs = None
		Xc = None
		Axis = None
		Theta = None
	elif key.char in "rR":
		print "Random samples"
		Xr = randomSample(Xg,5)
		print Xr
	elif key.char in "fF":
		print "Ellipse fitting"
		ellipseFit(Xg)
	elif key.char in 'aA':
		tkMessageBox.showinfo("iCircleFit", "Example of circle fitting by joohwi@cs.unc")	
	elif key.char in "qQ":
		exit(0)
			
def pointAdded():
	global Xg
	if (Xg is None):
		Xg = np.asarray([(mouse.x, mouse.y)])
	else:
		Xg = np.vstack((Xg, (mouse.x,mouse.y)))
		
def optimize():
	global Xg,Xr
	if (Xg is None or len(Xg) < 10):
		return
	K = 11
	for i in range(0,1000):
		Xr = randomSample(Xg,K)
		C = ellipseFit(Xr)
		T = 0.001
		nInliers = 0
		for x in Xg:
			cost = ellipseCostFunc(C,x)
			if (ellipseCostFunc(C,x) < T):
				nInliers += 1
		ratio = nInliers/float(len(Xg))
		if (ratio > 0.9):
			print "Best fit found"
			break

def mouseClicked():
	global Xg,Xs,Xi,Xr,X0
	pointAdded()

def test1(opts, argv):
	"Plotting test"
	Xs = op.fmin(f, np.asarray([0,0,1]), xtol=1e-10, full_output=0, disp=0)
	plotCircle(argv[0], Xi, Xs)


def main(opts, argv):
	run()

if __name__ == "__main__":
	parser = OptionParser(usage="%prog ")
	(opts, args) = parser.parse_args()

	main(opts, args)
