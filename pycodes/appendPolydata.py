#!/home/joohwi/python


import vtk,sys,niral as nu
from optparse import OptionParser

def main(opts, args):
  appender = vtk.vtkAppendPolyData()
  for f in args:
    o = nu.readVTK(f)
    appender.AddInput(o)
  appender.Update()
  nu.writeVTK(opts.output, appender.GetOutput())

if (__name__ == "__main__"):
  parser = OptionParser(usage="appendPolydata.py [options] input1 input2 ...")
  parser.add_option("-o", "--output", dest="output", help="output filename", default="output_appended.vtp")
  (opts, args) = parser.parse_args()
  if (len(args) == 0):
    parser.print_help()
  else:
    main(opts, args)

