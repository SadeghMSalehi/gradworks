#!/home/joohwi/python

from optparse import OptionParser
import vtk
import niral as nu
import numpy as np

def main(opts, args):
  inputTransform = np.array(opts.transform.split(","), dtype='double').reshape((4,4))
  transformMatrix = vtk.vtkMatrix4x4()
  for i in range(0,4):
    for j in range(0,4):
      transformMatrix.SetElement(i,j,inputTransform[i,j])

  print "InputTransform\n", inputTransform, "vtkMatrix\n", transformMatrix
  transform = vtk.vtkTransform()
  transform.SetMatrix(transformMatrix)
  transformFilter = vtk.vtkTransformPolyDataFilter()
  transformFilter.SetTransform(transform)
  for (inputId, input) in enumerate(args):
    obj = nu.readVTK(input)
    transformFilter.SetInput(obj)
    transformFilter.Update()
    objOut = transformFilter.GetOutput()
    if (opts.output.find("%") > 0):
      nu.writeVTK(opts.output % (inputId), objOut)
    else:
      nu.writeVTK(opts.output, objOut)
    
if __name__ == "__main__":
  parser = OptionParser(usage="transformPolydata.py [options] input1 input2 ...")
  parser.add_option("-o", "--output", dest="output", help="output transform polydata with index holder (%02d)", default="transformdOutput-%02d.vtp")
  parser.add_option("-t", "--transform", dest="transform", help="a comma separated array representing 4x4 matrix with row-major")
  (opts, args) = parser.parse_args()

  if (len(args) == 0):
    parser.print_help()
  else:
    main(opts,args)
