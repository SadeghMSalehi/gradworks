#!/home/joohwi/python

import sys
from vtkIOPython import *
from vtkFilteringPython import *
from vtkGraphicsPython import *
from vtkCommonPython import *
import csv
import niral as nu
from optparse import OptionParser
import numpy as np
import mdp
import math
import json

def loadJSON(f):
  print "Loading JSON:", f
  return json.loads(open(f).read())

def createCube(vertexList):
  cs = vtkCubeSource()
  cs.SetXLength(10)
  cs.SetYLength(10)
  cs.SetZLength(10)
  cs.Update()
  cp = vtkCleanPolyData()
  cp.SetInput(cs.GetOutput())
  cp.PointMergingOn()
  cp.ConvertStripsToPolysOn()
  cp.ConvertPolysToLinesOn()
  cp.ConvertLinesToPointsOn()
  cp.Update()
  cb = cp.GetOutput()
  vpts = cb.GetPoints()
  idx = [ 0,1,3,2,4,7,5,6 ]
  for (j,vertex) in enumerate(vertexList):
    vpts.SetPoint(idx[j],vertex[0],vertex[1],vertex[2])
  cb.SetPoints(vpts)
  return cb

# This code create bounding box as a geometry object
# Currently, this generates a cube
#
def buildBoundingBox(fname,fout):
  obbList = loadJSON(fname)
  appender = vtkAppendPolyData()
  for (id,label) in enumerate(obbList):
    vertexList = [];
    for vertex in sorted(obbList[label]['obb_vertices']):
      obbVertices = obbList[label]['obb_vertices'][vertex]
      vertexList.append(obbVertices)
    cubePolyData = createCube(vertexList)
    nu.writeVTK(fout % (id+1), cubePolyData) 
    appender.AddInput(cubePolyData)
  appender.Update()
  nu.writeVTK(fout%(0),appender.GetOutput())
      
# createEllipse takes center as [x,y,z], 
#   axesDirection as eigenvector columns, 
#   and axesLength as [rx,ry,rz]
def createEllipse(center,axesDirection,axesLength):
  scale = axesLength
  tx = vtkTransform()
  tx.PostMultiply()
  tx.Scale(scale[0],scale[1],scale[2])
  rotMat = vtkMatrix4x4()
  for i in range(0,3):
    for j in range(0,3):
      rotMat.SetElement(j,i,obbInfo['obb_rotation'][i*3+j])
  rotMat.SetElement(3,3,1)
  # skip rotation for a while
  tx.Concatenate(rotMat)
  translate = center
  logger.info("center: %s" % (translate))
  #translate = obbInfo['obb_vertices']
  tx.Translate(translate[0],translate[1],translate[2]);
  txf = vtkTransformPolyDataFilter()
  txf.SetTransform(tx)
  sphereSource = vtkSphereSource()
  sphereSource.SetRadius(1);
  sphereSource.SetThetaResolution(16);
  sphereSource.SetPhiResolution(16);
  sphereSource.Update()
  sphere = sphereSource.GetOutput()
  txf.SetInput(sphere)
  txf.Update()

def orientedBoundingBoxToTransform(obbInfo, components):
  tx = vtkTransform()
  tx.PostMultiply()
  scale = map(lambda x: x/2.0,obbInfo['obb_length'])
  tx.Scale(scale[0],scale[1],scale[2])
  rotMat = vtkMatrix4x4()
  for i in range(0,3):
    for j in range(0,3):
      rotMat.SetElement(j,i,obbInfo['obb_rotation'][i*3+j])
  rotMat.SetElement(3,3,1)
  # skip rotation for a while
  tx.Concatenate(rotMat)
  vertexList = obbInfo['obb_vertices'].values()
  center = [0,0,0]
  for vertex in vertexList:
    for j in range(0,3):
      center[j] = center[j] + vertex[j]
  for j in range(0,3):
    center[j] = center[j]/8.0
  translate = center
  #translate = obbInfo['obb_vertices']
  tx.Translate(translate[0],translate[1],translate[2]);
  components.append(scale)
  components.append(obbInfo['obb_rotation'])
  components.append(translate)
  return tx


def vtkMatrixToArray(mat):
  ''' convert vtk4x4 matrix into python array'''
  rowArray = []
  for i in range(0,4):
    for j in range(0,4):
      rowArray.append(mat.GetElement(i,j))
  return rowArray

def measureAngleBetweenAxis(axisDict):
  for key in axisDict.keys():
    for s in range(0,3):
      axisList = axisDict[key]
      axisList = axisList[s:len(axisList):3]
      N = len(axisList)
      M = np.zeros((N,N))
      for (i,axisMat1) in enumerate(axisList):
        for (j,axisMat2) in enumerate(axisList):
          if (i == j): 
            continue
          dot = abs(axisMat1.transpose().dot(axisMat2))
          dot = 1 if dot > 1 else 0 if dot < -1 else dot
          angle = math.degrees(math.acos(dot))
          M[i,j] = angle
      np.savetxt("AngleMeasure_" + key + "_%d.txt" % (s), M, fmt="%5.3f")

def buildParameters(jsonInputs):
  paramSet = []
  axisDict = {}
  for input in jsonInputs:
    obbList = loadJSON(input)
    allParam = []
    for (id,label) in enumerate(obbList):
      obbInfo = obbList[label]
      axisMatrix = np.matrix(obbInfo['eigenvectors']).reshape([3,3])
      print obbInfo['eigenvectors']
      print obbInfo['obb_rotation']
      if (label not in axisDict.keys()):
        axisDict[label] = [ ]
      axisDict[label].append(axisMatrix[:,0])
      axisDict[label].append(axisMatrix[:,1])
      axisDict[label].append(axisMatrix[:,2])
      tx = orientedBoundingBoxToTransform(obbInfo, [])
      txMat = tx.GetMatrix()
      txArr = vtkMatrixToArray(tx.GetMatrix())
      #print txArr[0:12], txMat
      allParam.extend(txArr[0:12])
    paramSet.append(allParam)
  paramMatrix = np.matrix(paramSet)
  np.savetxt("TransformParameter.txt", paramMatrix)
  measureAngleBetweenAxis(axisDict)

def convertFormat(opts, args):
  allTransforms = []
  for (jsonId, jsonfile) in enumerate(args):
    obbList = loadJSON(jsonfile)
    partsTransforms = []
    for (id,label) in enumerate(obbList):
      obbInfo = obbList[label]
      transforms = []
      tx = orientedBoundingBoxToTransform(obbInfo, transforms)
      partsTransforms.append(transforms)
    allTransforms.append(partsTransforms)
    #outfile = open(opts.output % jsonId, "w")
    #outfile.write(json.dumps(partsTransforms))
    #outfile.close()
  outfile = open(opts.output, "w")
  outfile.write(json.dumps(allTransforms))
  outfile.close()

def alignModels(opts, args):
  alltransforms = loadJSON(args[0])
  newTransforms = []
  for subj in alltransforms:
    subjTransforms = []
    T0 = np.array(subj[1][2])
    R0t = np.matrix(subj[1][1]).reshape(3,3).transpose()
    for part in subj:
      R = np.matrix(part[1]).reshape(3,3)
      Rx = R
      T = np.array(part[2])
      Tx = T - T0
      part[1] = Rx.ravel().tolist()[0]
      part[2] = Tx.tolist()
      subjTransforms.append([part[0], part[1], part[2]])
    newTransforms.append(subjTransforms)
  outfile = open("allTransformsAligned.txt", "w")
  outfile.write(json.dumps(newTransforms))
  outfile.close()

def createTransform(S, R, T):
  tx = vtkTransform()
  tx.PostMultiply()
  tx.Scale(S[0], S[1], S[2])
  rotMat = vtkMatrix4x4()
  for i in range(0,3):
    for j in range(0,3):
      rotMat.SetElement(j,i,R[j,i])
  rotMat.SetElement(3,3,1)
  # skip rotation for a while
  tx.Concatenate(rotMat)
  tx.Translate(T[0],T[1],T[2])
  return tx


def buildEllipses2(opts, args):
  alltransforms = loadJSON(args[0])
  for (subjId, subj) in enumerate(alltransforms):
    print "Processing", subjId
    appender = vtkAppendPolyData()
    for part in subj:
      try:
        S = part[0]
        R = np.matrix(part[1]).reshape(3,3).transpose()
        T = part[2]
      except:
        print "Error:", part
        return
      print S, R, T
      txf = vtkTransformPolyDataFilter()
      txf.SetTransform(createTransform(S, R, T))
      sphereSource = vtkSphereSource()
      sphereSource.SetRadius(1);
      sphereSource.SetThetaResolution(16);
      sphereSource.SetPhiResolution(16);
      sphereSource.Update()
      sphere = sphereSource.GetOutput()
      txf.SetInput(sphere)
      txf.Update()
      #nu.writeVTK(opts.output%(id+1), txf.GetOutput()) 
      appender.AddInput(txf.GetOutput())
    appender.Update()
    fileOut = opts.output % (subjId)
    print "Saving", fileOut
    nu.writeVTK(opts.output % (subjId), appender.GetOutput())


def runPCA(paramFile):
  data = np.loadtxt(paramFile)
  y = mdp.pca(data, reduce=True)
  print y

def viewGeometry(opts, args):
  geom = loadJSON(args[0])
  for (subjId, subj) in enumerate(geom):
    print "Subj:", subjId
    for parts in subj:
      print parts[0], parts[2]

def main():
  parser = OptionParser(usage="%prog [options] bounding_box_json.txt")
  parser.add_option("-b", dest="bboxOut", help="specify a filename for bounding box output with id holder (%d)", default="")
  parser.add_option("-E", dest="ellipseOut2", help="generate an ellipse models from the first arg and save into the file given by -o", action="store_true")
  parser.add_option("-p", dest="paramOut", help="specify a filename for parameter output with id holder (%d)", default="")
  parser.add_option("-a", dest="analysis", help="run PCA analysis", default=False, action="store_true")
  parser.add_option("-c", dest="convert", help="convert json to transformation", action="store_true")
  parser.add_option("-o", dest="output", help="output with place holder(%s)")
  parser.add_option("-l", dest="alignment", help="align to center ellipse", action="store_true")
  parser.add_option("-v", dest="view", help="view part geometry information", action="store_true")

  (opts,args) = parser.parse_args()
  if (len(args) == 0):
    parser.print_help()
  else:
    print "Running ..."
    if (opts.paramOut != ""):
      buildParameters(args)
    elif (opts.bboxOut != ""):
      buildBoundingBox(args[0], opts.bboxOut)
    elif (opts.convert):
      convertFormat(opts, args)
    elif (opts.alignment):
      alignModels(opts, args) 
    if (opts.ellipseOut2):
      print "Build Ellipses ..."
      buildEllipses2(opts, args)
    if (opts.analysis):
      runPCA(args[0])
    if (opts.view):
      viewGeometry(opts, args)

if (__name__ == "__main__"):
  main()


