#!/tools/ParaView-Build/bin/pvpython

from vtk import *
from optparse import OptionParser
import math

import numpy as np
import scipy as sp
import scipy.optimize as op
import scipy.ndimage.filters as spnf
import scipy.sparse as scsp
import networkx as nx
import glob
import matplotlib.pyplot as ppl

def readVTP(f):
	if (f.endswith(".vtp")):
		r = vtkXMLPolyDataReader()
		r.SetFileName(f)
		r.Update()
		return r.GetOutput()
	elif (f.endswith(".vtk")):
		r = vtkPolyDataReader()
		r.SetFileName(f)
		r.Update()
		return r.GetOutput()
	elif (f.endswith(".vtu")):
		r = vtkXMLUnstructuredGridReader()
		r.SetFileName(f)
		r.Update()
		return r.GetOutput()

def writeVTK(f, p):
	if (f.endswith("vtu")):
		w = vtkXMLUnstructuredGridWriter()
		w.SetFileName(f)
		w.SetInput(p)
		w.Update()
	elif (f.endswith("vtp")):
		w = vtkXMLPolyDataWriter()
		w.SetFileName(f)
		w.SetInput(p)
		w.Update()
	else:
		w = vtkPolyDataWriter()
		w.SetFileName(f)
		w.SetInput(p)
		w.Update()

def enumPoints(l,s):
	pMat = np.zeros([l.GetNumberOfIds(),3])
	for x in range(0, l.GetNumberOfIds()):
		i = l.GetId(x)
		p = s.GetPoint(i)
		pMat[x,0] = p[0]
		pMat[x,1] = p[1]
		pMat[x,2] = p[2]
	return pMat

def enumPoints2(pset, ids):
	p = vtkPoints()
	for i in range(0, ids.GetNumberOfIds()):
		p.InsertNextPoint(pset.GetPoint(ids.GetId(i)))
	return p

def updatePoints(l,s,pIn):
	for x in range(0, l.GetNumberOfIds()):
		i = l.GetId(x)

		pUpdate = [0,0,0]
		pUpdate[0] = pIn[x,0]
		pUpdate[1] = pIn[x,1]
		pUpdate[2] = pIn[x,2]
		s.SetPoint(i, pUpdate)

def enumLines(p, func):
	iterLine = p.GetLines()
	iterLine.InitTraversal()
	idList = vtkIdList()
	lineId = 1

	while (lineId > 0):
		lineId = iterLine.GetNextCell(idList)
		aLine = enumPoints(idList, p.GetPoints())

		if (len(aLine) > 100):
			func(aLine, idList, lineId, p)

def writeLineAsMAT(opts, l):
	s = l[opts.id]
	f = open(opts.output, "w")
	for p in s:
		f.write("%f\t%f\t%f\n" % (p[0], p[1], p[2]))
	f.close()

def optimizeSphere3d(data):
	"Fit a 3d sphere into data, initial is given in x0 (x,y,z,r)"
	def f(x):
		return np.power(np.linalg.norm(data-x[0:2].sum())-x[3],2)
	c = data.mean(0)	
	r = (data.max(0)-data.min(0)).min()
	x0 = (c[0],c[1],c[2],r)
	print x0
	xs = op.fmin(f,x0,xtol=1e-10,disp=0)
	return xs

def sphereAt(xs):
	s = vtk.vtkSphereSource()
	s.SetCenter(xs[0],xs[1],xs[2])
	s.SetRadius(xs[3])
	s.Update()
	return s.GetOutput()

def fitSpheres(opts,argv):
	"Compute fitting sphere for vtkPolyData with lines"
	p = readVTP(argv[0])
	a = vtk.vtkAppendPolyData()
	for l in lineGenerator(p.GetLines()):
		pSet = enumPoints(l,p.GetPoints())
		"return (x,y,z,r)"
		xs = optimizeSphere3d(pSet)
		a.AddInput(sphereAt(xs))
	a.Update()
	writeVTK(argv[1], a.GetOutput())


# lps: line-points
# ipm: id - point mapping
# lid: line id
# poly: polygonal mesh
def computeCurvature(lps, ipm, lid, poly):
	kArray = poly.GetPointData().GetScalars()

	r_lps = np.roll(lps, 1, 0)
	d1 = r_lps - lps
	r_d1 = np.roll(d1, 1, 0)
	d2 = r_d1 - d1

	for i in range(0, ipm.GetNumberOfIds()):
		dp = d1[i,:]
		d2p = d2[i,:]

		A = dp[1]*d2p[2] - dp[2]*d2p[1]
		B = dp[2]*d2p[0] - dp[1]*d2p[2]
		C = dp[0]*d2p[1] - dp[1]*d2p[0]
		
		denom = math.pow(dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2], 1.5)
		k = math.sqrt(A*A + B*B + C*C) / denom

		kArray.SetValue(ipm.GetId(i), k)
	return


class OBB:
	corner = np.zeros(3)
	m1 = np.zeros(3)
	m2 = np.zeros(3)
	m3 = np.zeros(3)
	sz = np.zeros(3)

	def __repr__(self):
		return "(%f,%f,%f) (%f,%f,%f) (%f,%f,%f) (%f,%f,%f) (%f,%f,%f)" % (self.corner[0], self.corner[1], self.corner[2], self.m1[0],self.m1[1],self.m1[2],self.m2[0],self.m2[1],self.m2[2],self.m3[0],self.m3[1],self.m3[2],self.sz[0],self.sz[1],self.sz[2])

	def MakeCube(self, pset):
		iStart = pset.GetNumberOfPoints()
		pset.InsertNextPoint(self.corner)
		pset.InsertNextPoint(self.corner+self.m1)
		pset.InsertNextPoint(self.corner+self.m2)
		pset.InsertNextPoint(self.corner+self.m1+self.m2)
		pset.InsertNextPoint(self.corner+self.m3)
		pset.InsertNextPoint(self.corner+self.m3+self.m1)
		pset.InsertNextPoint(self.corner+self.m3+self.m2)
		pset.InsertNextPoint(self.corner+self.m3+self.m1+self.m2)
		ids = vtkIdList()
		for i in range(iStart, iStart + 8):
			ids.InsertNextId(i)
		return ids

def lineToCubes(poly):
	idList = vtkIdList()

	iterLine = poly.GetLines()
	iterLine.InitTraversal()
	keepIter = 1

	o = []

	pst = vtk.vtkPoints()

	g = vtk.vtkUnstructuredGrid()
	g.SetPoints(pst)

	while (keepIter > 0):
		keepIter = iterLine.GetNextCell(idList)
		linePts = enumPoints2(poly.GetPoints(), idList)

		obb = OBB()
		vtk.vtkOBBTree.ComputeOBB(linePts, obb.corner, obb.m1, obb.m2, obb.m3, obb.sz)

		cube = obb.MakeCube(pst)
		g.InsertNextCell(vtk.VTK_VOXEL, cube)


		#p = vtkPolyData()
		#p.SetPoints(poly.GetPoints())
		
		#c = vtkCellArray()
		#c.InsertNextCell(idList)

		#p.SetLines(c)

		#p.ComputeBounds()

		#bb = vtkBoundingBox()
		#bb.SetBounds(p.GetBounds())

	return g

def conv2mat(poly, fout):
	iterLine = poly.GetLines()	
	iterLine.InitTraversal()

	keepIter = 1
	idl = vtkIdList()

	f = open(fout, "w")

	keepIter = iterLine.GetNextCell(idl)
	while (keepIter > 0):
		pts = poly.GetPoints()

		print "Writing # of points: %d" % (idl.GetNumberOfIds())
		for i in range(0, idl.GetNumberOfIds()):
			p = pts.GetPoint(idl.GetId(i))
			f.write("%f %f %f\n" % (p[0], p[1], p[2]))

		keepIter = iterLine.GetNextCell(idl)


	f.close()

def smoothLine(lps):
	H = np.asarray([ 0.0545, 0.2442, 0.4026, 0.2442, 0.0545 ])

	outlps = np.zeros([lps.shape[0]-4,3])
	for d in range(0,3):
		outlps[:,d] = np.convolve(lps[:,d], H, 'valid')
	
	return outlps

def firstLineAsArray(polyLine, ids):
	iterLine = polyLine.GetLines()	
	iterLine.InitTraversal()

	keepIter = 1
	idl = vtkIdList()

	keepIter = iterLine.GetNextCell(idl)

	x = np.zeros((idl.GetNumberOfIds(), 3))
	for i in range(0, idl.GetNumberOfIds()):
		pId = idl.GetId(i)
		p = polyLine.GetPoints().GetPoint(pId)
		x[i,:] = p
		ids.append(pId)
	return x

def computeCurvatureForLine(lps):
	if (lps.shape[0] < 5):
		return None

	r_lps = np.roll(lps, 1, 0)
	d1 = r_lps - lps
	r_d1 = np.roll(d1, 1, 0)
	d2 = r_d1 - d1

	kout = np.zeros(lps.shape[0])
	for i in range(0, lps.shape[0]):
		dp = d1[i,:]
		d2p = d2[i,:]

		A = dp[1]*d2p[2] - dp[2]*d2p[1]
		B = dp[2]*d2p[0] - dp[1]*d2p[2]
		C = dp[0]*d2p[1] - dp[1]*d2p[0]
		
		denom = math.pow(dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2], 1.5)
		k = math.sqrt(A*A + B*B + C*C) / denom
		kout[i] = k

	return kout[2:len(kout)-2]

def bracketPeaks(kLine, firstTh, lastTh):
	firstPeakId = 0
	lastPeakId  = len(kLine) - 1 

	dKLine = kLine - np.roll(kLine, 1, 0)
	for (idx, dK) in enumerate(dKLine):
		if (dK < firstTh):
			firstPeakId = idx + 2
			break

	for (idx, dK) in enumerate(dKLine[::-1]):
		if (dK > lastTh):
			lastPeakId = len(kLine) - idx - 2
			break

	return (firstPeakId, lastPeakId, dKLine)

	
def filterStreams(poly):
	newLines = vtkCellArray()
	for idl in lineGenerator(poly.GetLines()):
		x = np.zeros((idl.GetNumberOfIds(), 3))
		if (idl.GetNumberOfIds() < 80):
#				print "# of ID %d" % (idl.GetNumberOfIds())
			continue

		for i in range(0, idl.GetNumberOfIds()):
			pId = idl.GetId(i)
			p = poly.GetPoints().GetPoint(pId)
			x[i,:] = p

		kLine = computeCurvatureForLine(x)
		if (kLine.max() < 1000):
#				print "min of K = %f" % (kLine.max())
			continue
	
		(firstPeakId, lastPeakId, dKLine) = bracketPeaks(kLine, -400, 400)

		if (firstPeakId > lastPeakId or (lastPeakId - firstPeakId) < 80):
#				print "bad bracket[%d:%d]" % (firstPeakId, lastPeakId)
			continue

		newLineIds = vtk.vtkIdList()
		print "Processing [%d:%d]" % (firstPeakId, lastPeakId)
		for idx in range(firstPeakId, lastPeakId):
			newLineIds.InsertNextId(idl.GetId(idx))

		print kLine.max()
		newLines.InsertNextCell(newLineIds)

	pOut = vtk.vtkPolyData()
	pOut.SetPoints(poly.GetPoints())
	pOut.SetLines(newLines)

	print "# of lines: %d" % (newLines.GetNumberOfCells())
	return pOut

def lineGenerator(polyLines):
	idl = vtkIdList()
	polyLines.InitTraversal()
	keepIter = polyLines.GetNextCell(idl)
	while (keepIter > 0):
		try:
			yield idl
		finally:
			keepIter = polyLines.GetNextCell(idl)
			
def idListToPoints(idList, points):
	pMat = np.zeros([idList.GetNumberOfIds(),3])
	for x in range(0, idList.GetNumberOfIds()):
		i = idList.GetId(x)
		p = points.GetPoint(i)
		pMat[x,0] = p[0]
		pMat[x,1] = p[1]
		pMat[x,2] = p[2]
	return pMat

class AABB:
	"Axis-aligned bounding box"
	P = np.asarray([np.inf, np.inf, np.inf])
	Q = np.asarray([-np.inf, -np.inf, -np.inf])
	E = [0,0,0]
	V = 0.0

	def __init__(self, *args, **kwargs):
		if (len(args) == 1):
			for p in args[0]:
				self.AddPoint(p)
		elif (len(args) == 2):
			P = args[0]
			Q = args[1]
			self.P = np.minimum(P, Q)
			self.Q = np.maximum(P, Q)
			self.E = Q-P
			self.V = self.E[0]*self.E[1]*self.E[2]

	def AddPoint(self, pts):
		self.P = np.minimum(self.P, pts)
		self.Q = np.maximum(self.Q, pts)
		self.E = self.Q-self.P
		self.V = self.E[0] * self.E[1] * self.E[2]

	def Intersect(self, o):
		"AABB-AABB intersection. Return its volume"
		newE = np.minimum(self.Q, o.Q) - np.maximum(self.P, o.P)
		if (newE[0] <= 0 or newE[1] <= 0 or newE[2] <= 0):
			return 0
		else:
			v = newE[0]*newE[1]*newE[2]
			minV = min(self.V, o.V)
			v = v / minV if minV > 0 else 0
			return v

	def __repr__(self):
		return "(%f,%f,%f) - (%f,%f,%f)" % \
			(self.P[0], self.P[1], self.P[2], \
			self.Q[0], self.Q[1], self.Q[2])

	def toArray(self):
		return np.concatenate([self.P, self.Q])

def addCellIntArray(poly, name, arr):
	v = vtkIntArray()
	v.SetName(name)
	v.SetNumberOfValues(len(arr))
	for (i,value) in enumerate(arr):
		v.SetValue(i,value)
	poly.GetCellData().AddArray(v)

def glyphAABB(aabbs):
	pset = vtk.vtkPoints()
	g = vtk.vtkUnstructuredGrid()
	g.SetPoints(pset)
	cellIds = vtkIntArray()
	cellIds.SetName("CellIds")
	for (i,aabb) in enumerate(aabbs):
		ids = vtkIdList()
		iStart = pset.GetNumberOfPoints()
		pset.InsertNextPoint(aabb.P)
		pset.InsertNextPoint(aabb.P+(aabb.E[0],0,0))
		pset.InsertNextPoint(aabb.P+(0,aabb.E[1],0))
		pset.InsertNextPoint(aabb.P+(aabb.E[0],aabb.E[1],0))
		pset.InsertNextPoint(aabb.P+(0,0,aabb.E[2]))
		pset.InsertNextPoint(aabb.P+(aabb.E[0],0,aabb.E[2]))
		pset.InsertNextPoint(aabb.P+(0,aabb.E[1],aabb.E[2]))
		pset.InsertNextPoint(aabb.P+(aabb.E[0],aabb.E[1],aabb.E[2]))
		ids = vtkIdList()
		for i in range(iStart, iStart + 8):
			ids.InsertNextId(i)
		g.InsertNextCell(vtk.VTK_VOXEL, ids)
		cellIds.InsertNextValue(i)
	g.GetCellData().AddArray(cellIds)
	return g

def clusterByAABB(poly,opts,argv):
	aabbs = []
	aabbsOut = np.zeros((poly.GetLines().GetNumberOfCells(),6))
	for (i,idl) in enumerate(lineGenerator(poly.GetLines())):
		pts = idListToPoints(idl, poly.GetPoints())
		aabb = AABB(pts)
		aabbs.append(aabb)
		aabbsOut[i,:] = aabb.toArray()
	overlaps = np.zeros((len(aabbs), len(aabbs)))
	for (i1, ab1) in enumerate(aabbs):
		for i2 in range(i1, len(aabbs)):
			if (i1 == i2):
				overlaps[i1,i2] = 1
			else:
				ab2 = aabbs[i2]
				overlaps[i1,i2] = ab1.Intersect(ab2)
				overlaps[i2,i1] = overlaps[i1,i2]
	threshold = np.vectorize(lambda x: 1 if x > opts.overlap else 0, otypes=[np.int32])
	adjacency = threshold(overlaps)
	conComps = scsp.cs_graph_components(adjacency)
	addCellIntArray(poly, "VortexCluster", conComps[1])
	aabbGlyphs = glyphAABB(aabbs)
	addCellIntArray(aabbGlyphs, "VortexCluster", conComps[1])
	if (opts.aabbOut != ""):
		np.savetxt(opts.aabbOut, aabbsOut)
	if (opts.aabbLabels != ""):
		np.savetxt(opts.aabbLabels, conComps[1])
	return (poly, aabbGlyphs)

def createOBB(poly):
	obb = vtk.vtkOBBTree()
	obb.SetMaxLevel(10)
	obb.SetNumberOfCellsPerNode(5)
	obb.AutomaticOn()
	sfilt = vtk.vtkSpatialRepresentationFilter()
	sfilt.SetSpatialRepresentation(obb)
	sfilt.SetInput(poly)
	sfilt.Update()
	return sfilt.GetOutput()

def trackAABB(opts, argv):
	"AABB Tracking"

	""" In case there's no label file
	flist = glob.glob(argv[0])
	for f in flist:
		p = readVTP(f)
		k = vtk.mutable(0)
		a = p.GetCellData().GetArray("VortexCluster", k)
		x = np.zeros(a.GetNumberOfTuples())
		for i in range(0, a.GetNumberOfTuples()):
			v = a.GetValue(i)
			x[i] = v
		np.savetxt(f.replace("AABB.vtu", "Cluster_labels.txt"), x, fmt="%d")
	"""
	"Load ABB"
	g = nx.Graph()
	flist = sorted(glob.glob(argv[1]))
	labels = []
	for (t,f) in enumerate(flist):
		L = np.genfromtxt(f, dtype=int)
		for i in range(0, L.max()+1):
			g.add_node("%d.%d" % (t,i), time=t)
		labels.append(L)
		if (t >= opts.maxT):
			break
	flist = sorted(glob.glob(argv[0]))
	prevaabbs = []
	for (t,f) in enumerate(flist):
		print "Processing time #%d" % (t)
		aabb = np.genfromtxt(f)
		aabbs = []
		for ab in aabb:
			aabbObj = AABB(ab[0:3],ab[3:6])
			aabbs.append(aabbObj)
		if (t > 0):
			for (i,x) in enumerate(prevaabbs):
				ci = "%d.%d" % (t-1,labels[t-1][i])
				for (j,y) in enumerate(aabbs):
					cj = "%d.%d" % (t,labels[t][j])
					v = x.Intersect(y)
					if (v > opts.overlap):
						g.add_edge(ci,cj)
		prevaabbs = aabbs
		if (t >= opts.maxT):
			break
	nx.write_adjlist(g, argv[2])
	if (opts.subgraph != ""):
		C = nx.connected_component_subgraphs(g)
		for i,c in enumerate(C):
			nx.write_adjlist(c, opts.subgraph % (i))

def timeflow(opts, argv):
	"""
	Read cluster tracking results and aggregate into a single file
	"""
	g = nx.read_adjlist(argv[0])
	f = sorted(glob.glob(opts.aabbIn))
	N = map(lambda x: map(int,x.split(".")), nx.nodes(g))
	C = dict((t,set()) for t in map(lambda x:x[0], N))
	for (t,l) in N: C[t].add(l)
	newMesh = vtk.vtkPolyData()
	newLines = vtk.vtkCellArray()
	newPoints = vtk.vtkPoints()
	newTimeData = vtkIntArray()
	newTimeData.SetName("TimeStep")
	for t in C:
		p = readVTP(f[t])
		"Filter cluster labels"
		a = p.GetCellData().GetArray("VortexCluster", vtk.mutable(0))
		for (i,l) in enumerate(lineGenerator(p.GetLines())):
			c = a.GetValue(i)
			if c in C[t]:
				newLine = vtkIdList()
				for j in range(0, l.GetNumberOfIds()):
					newLine.InsertNextId(newPoints.GetNumberOfPoints())
					newPoints.InsertNextPoint(p.GetPoint(l.GetId(j)))
				newLines.InsertNextCell(newLine)
				newTimeData.InsertNextValue(t)
	newMesh.GetCellData().SetScalars(newTimeData)
	newMesh.SetPoints(newPoints)
	newMesh.SetLines(newLines)
	writeVTK(opts.output, newMesh)
			
def main(opts, argv):
	if (opts.curvature):
		p = readVTP(argv[0])
		k = vtkDoubleArray()
		k.SetNumberOfComponents(1)
		k.SetNumberOfValues(p.GetNumberOfPoints())
		k.SetName("Line Curvature")
		p.GetPointData().SetScalars(k)
		enumLines(p, computeCurvature)
		writeVTK(argv[1], p)
	elif (opts.cube):
		p = readVTP(argv[0])
		o = lineToCubes(p)
		fout = argv[1]
		writeVTK(fout, o)
	elif (opts.mat):
		p = readVTP(argv[0])
		conv2mat(p, argv[1])	
	elif (opts.aabbTracking):
		trackAABB(opts, argv)
	elif (opts.timeflow):
		timeflow(opts,argv)
	elif (opts.fitspheres):
		fitSpheres(opts,argv)
	elif (opts.aabb):
		p = readVTP(argv[0])
		fout = argv[1]
		(pOut, pAABB) = clusterByAABB(p, opts, argv)
		writeVTK(fout, pOut)
		if (len(argv) > 2):
			writeVTK(argv[2], pAABB)
	elif (opts.obb):
		p = readVTP(argv[0])
		pOut = createOBB(p)
		writeVTK(argv[1], pOut)
	else:
		p = readVTP(argv[0])
		pOut = filterStreams(p)
		fout = argv[1]
		writeVTK(fout, pOut)

if (__name__ == "__main__"):
	parser = OptionParser(usage="%prog [options] vtp-file")
	parser.add_option("-i", dest="id", help="line id", type=int, default=0)
	parser.add_option("-o", dest="output", help="output file name", default="aline.mat")
	parser.add_option("-k", dest="curvature", help="compute line curvature values", action="store_true", default=False)
	parser.add_option("-c", dest="cube", help="compute bounding boxes", action="store_true", default=False) 
	parser.add_option("-m", dest="mat", help="conversion vtk2mat", action="store_true", default=False) 
	parser.add_option("--aabb", dest="aabb", help="clustering with AABB", action="store_true", default=False)
	parser.add_option("--overlap", dest="overlap", help="Overlapping threshold", default=0.8)
	parser.add_option("--aabbIn", dest="aabbIn", help="AABB point coordinates input", default="")
	parser.add_option("--aabbOut", dest="aabbOut", help="AABB point coordinates", default="")
	parser.add_option("--aabbLabels", dest="aabbLabels", help="AABB label per streamline", default="")
	parser.add_option("--obb", dest="obb", help="create obb tree using VTK", action="store_true", default=False)
	parser.add_option("--aabbtracking", dest="aabbTracking", help="AABB tracking", action="store_true", default=False)
	parser.add_option("--maxT", dest="maxT", help="maximum timestep for AABB tracking", type=int, default=100)
	parser.add_option("--subgraph", dest="subgraph", help="subgraph output of the tracking", default="")
	parser.add_option("--timeflow", dest="timeflow", help="[--aabbIn streamlines with VortexCluster attribute -o output.vtp graphComponent.txt] processing tracking result to show flow over time", action="store_true", default=False)
	parser.add_option("--fitspheres", dest="fitspheres", help="fit spheres for each streamline", action="store_true", default=False)

	(opts, argv) = parser.parse_args()
	main(opts, argv)
