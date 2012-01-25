#!/usr/bin/python

import vtk
import niralvtk as nv
import numpy as np
from optparse import OptionParser

class OBB:
	def __init__(self):
		self.corner = np.zeros(3)
		self.m1 = np.zeros(3)
		self.m2 = np.zeros(3)
		self.m3 = np.zeros(3)
		self.sz = np.zeros(3)
		self.center = np.zeros(3)

	def Center(self):
		return self.corner+(self.m1+self.m2+self.m3)/2.0
		
	def __repr__(self):
		return "(%f,%f,%f) (%f,%f,%f) (%f,%f,%f) (%f,%f,%f) (%f,%f,%f)" % (self.corner[0], self.corner[1], self.corner[2], self.m1[0],self.m1[1],self.m1[2],self.m2[0],self.m2[1],self.m2[2],self.m3[0],self.m3[1],self.m3[2],self.sz[0],self.sz[1],self.sz[2])

	def AxisAlign(self):
		(l1,l2,l3) = (np.linalg.norm(self.m1),np.linalg.norm(self.m2),np.linalg.norm(self.m3))
		(e1,e2,e3) = (self.m1/l1,self.m2/l2,self.m3/l3)
		(s1,s2,s3) = (np.dot((1,0,0),e1),np.dot((0,1,0),e2),np.dot((0,0,1),e3))
		(e1,e2,e3) = (s1*e1,s2*e2,s3*e3)
		r1 = RotateByAxis(e1,[1,0,0])
		r2 = RotateByAxis(e2,[0,1,0])
		return np.dot(r2,r1)

	def MakeCube(self, g):
		pset = g.GetPoints()
		iStart = pset.GetNumberOfPoints()
		corner = self.corner
		corner = corner-(self.corner+(self.m1+self.m2+self.m3)/2)
		pset.InsertNextPoint(corner)
		pset.InsertNextPoint(corner+self.m1)
		pset.InsertNextPoint(corner+self.m2)
		pset.InsertNextPoint(corner+self.m1+self.m2)
		pset.InsertNextPoint(corner+self.m3)
		pset.InsertNextPoint(corner+self.m3+self.m1)
		pset.InsertNextPoint(corner+self.m3+self.m2)
		pset.InsertNextPoint(corner+self.m3+self.m1+self.m2)
		ids = vtk.vtkIdList()
		for i in range(iStart, iStart + 8):
			ids.InsertNextId(i)
		g.InsertNextCell(vtk.VTK_VOXEL,ids)

def CellIterator(cells):
  idl = vtkIdList()
  cells.InitTraversal()
  keepIter = cells.GetNextCell(idl)
  while (keepIter > 0): 
    try:
      yield idl 
    finally:
      keepIter = cells.GetNextCell(idl)

def RotateByAxis(Axis1,Axis2):
	CAxis = np.cross(Axis1, Axis2)
	if np.linalg.norm(CAxis) > 0:
		CAxis = CAxis / np.linalg.norm(CAxis)
	CTheta = np.arccos(np.dot(Axis1, Axis2))
	A1 = np.asarray([[0, -CAxis[2], CAxis[1]],[CAxis[2], 0, -CAxis[0]], [-CAxis[1], CAxis[0], 0]])
	Rotation1 = np.eye(3) + A1 * np.sin(CTheta)+ np.dot(A1,A1 * (1 - np.cos(CTheta)))
	return Rotation1

def vtkTransformMatrix(m):
	M = vtk.vtkMatrix4x4()
	for i in range(0,3):
		for j in range(0,3):
			M.SetElement(i,j,m[i][j])
	for i in range(0,3):
		M.SetElement(3,i,0)
		M.SetElement(i,3,0)
	M.SetElement(3,3,1)
	return M
	
def concomp(opts, argv):
	p = nv.readVTK(argv[0])
	c = vtk.vtkPolyDataConnectivityFilter()
	c.SetInput(p)
	c.SetExtractionModeToLargestRegion()
	c.Update()
	d = vtk.vtkCleanPolyData()
	d.SetInput(c.GetOutput())
	d.Update()
	p = d.GetOutput()
	nv.writeVTK(argv[1],p)

def transform(p, T):
	t = vtk.vtkTransform()
	t.SetMatrix(vtkTransformMatrix(T))
	tx = vtk.vtkTransformPolyDataFilter()
	tx.SetTransform(t)
	tx.SetInput(p)
	tx.Update()
	return tx.GetOutput()

def translateToMean(P,m,s):
	for i in range(0,P.GetNumberOfPoints()):
		p = P.GetPoint(i)
		q = [0,0,0]
		for j in range(0,3):
			q[j] = p[j] + s*m[j]
		P.GetPoints().SetPoint(i,q)

def principalAxes(opts, argv):
	p = nv.readVTK(argv[0])
	for i in range(0,p.GetNumberOfPoints()):
		s = p.GetPoint(i)
		s = (-s[0],-s[1],s[2])
		p.GetPoints().SetPoint(i,s)
	#g = vtk.vtkUnstructuredGrid()
	#g.SetPoints(vtk.vtkPoints())
	pOBB = OBB()
	vtk.vtkOBBTree.ComputeOBB(p.GetPoints(), pOBB.corner, pOBB.m1, pOBB.m2, pOBB.m3, pOBB.sz)
	#pOBB.MakeCube(g)
	translateToMean(p,pOBB.Center(),-1)
	pR = pOBB.AxisAlign()
	print pR
	pOut = transform(p,pR)
	nv.writeVTK(argv[1], pOut)
	writeITKTransform(argv[2], pR, np.asarray((255,255,79)) - pOBB.Center(), pOBB.Center())

def icp(opts,args):
	v1 = nv.readVTK(args[0])
	v2 = nv.readVTK(args[1])
	icp = vtk.vtkIterativeClosestPointTransform()
	icp.SetSource(v2)
	icp.SetTarget(v1)
	icp.StartByMatchingCentroidsOn()
	icp.SetMaximumNumberOfIterations(500)
	icp.SetMaximumNumberOfLandmarks(20000)
	tx = icp.GetLandmarkTransform()
	tx.SetModeToRigidBody()
	icp.Update()
	tx = icp.GetLandmarkTransform()
	txf = vtk.vtkTransformPolyDataFilter()
	txf.SetInput(v2)
	txf.SetTransform(tx)
	txf.Update()
	v3 = txf.GetOutput()
	nv.writeVTK(args[2], v3)
	return

def obb(opts,argv):
	poly = nv.readVTK(argv[0])
	obb = vtk.vtkOBBTree()
	obb.SetMaxLevel(10)
	obb.SetNumberOfCellsPerNode(5)
	obb.AutomaticOn()
	sfilt = vtk.vtkSpatialRepresentationFilter()
	sfilt.SetSpatialRepresentation(obb)
	sfilt.SetInput(poly)
	sfilt.Update()
	nv.writeVTK(argv[1],sfilt.GetOutput())

def writeITKTransform(fname,m3x3,t,o):
	"must give inverse matrix"
	m = np.linalg.inv(m3x3)
	header = "#Insight Transform File V1.0\n# Transform 0\n"
	typeheader = "Transform: AffineTransform_double_3_3\n"
	params = "Parameters: %f %f %f %f %f %f %f %f %f %f %f %f\n" % \
		(m[0,0],m[0,1],m[0,2],m[1,0],m[1,1],m[1,2],m[2,0],m[2,1],m[2,2],t[0],t[1],t[2])
	fixedParams = "FixedParameters: %f %f %f\n" % (o[0],o[1],o[2])
	f = open(fname, "w")
	f.write(header + typeheader + params + fixedParams)
	f.close()

def normalHistogram(opts,argv):
	poly = nv.readVTK(argv[0])
	normfilt = vtk.vtkPolyDataNormals()
	normfilt.ConsistencyOn()
	normfilt.AutoOrientNormalsOn()
	normfilt.ComputeCellNormalsOn()
	normfilt.SetInput(poly)
	normfilt.SetFeatureAngle(90)
	normfilt.Update()
	poly = normfilt.GetOutput()
	nv.writeVTK("nrom.vtk", poly)
	normals = poly.GetCellData().GetArray("Normals", vtk.mutable(0))
	normMap = vtk.vtkUnstructuredGrid()
	normMap.SetPoints(vtk.vtkPoints())
	normPts = normMap.GetPoints()
	for i in range(0, normals.GetNumberOfTuples()):
		norm = [0,0,0]
		normals.GetTuple(i, norm)
		normPts.InsertNextPoint(norm)
	nv.writeVTK(argv[1], normMap)
	
def main(opts, argv):
	if (opts.concomp):
		concomp(opts, argv)
	elif opts.principalAxes:
		principalAxes(opts,argv)
	elif opts.rotate:
		rotate(opts,argv)
	elif opts.command != "":
		funcCommand = eval(opts.command)
		if funcCommand is not None:
			funcCommand(opts,argv)
		else:
			print opts.command + " is not a command"

if __name__ == "__main__":
	parser = OptionParser(usage="%prog [options] arg1 args2")
	parser.add_option("-o", "--output", dest="output", help="output file", default="")
	parser.add_option("-c", "--concomp", dest="concomp", help="connected component", action="store_true", default=False)
	parser.add_option("-p", "--principalAxes", dest="principalAxes", help="Principal Axes", action="store_true", default=False)
	parser.add_option("-r", "--rotate", dest="rotate", help="Rotate", action="store_true", default=False)
	parser.add_option("-C", "--command", dest="command", help="Command", default="")

	(opts, argv) = parser.parse_args()
	main(opts, argv)
