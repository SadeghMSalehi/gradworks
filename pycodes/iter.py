#!/usr/bin/python

from optparse import OptionParser
import niral
import os,glob,os.path

def launch(opts, args):
  cmd = " ".join(args)
  dirs = glob.glob(opts.dir)
  dirs.sort()
  for d in dirs:
    files = niral.findFilesByPattern(d, opts.pattern)
    for f in files:
      ncmd = cmd.replace("@", os.path.basename(f))
      print niral.red(os.path.dirname(f) + "$ " + ncmd)
      if (not opts.istest):
        cwd = os.getcwd()
        os.chdir(os.path.dirname(f))
        os.system(ncmd)
        os.chdir(cwd)

if __name__ == "__main__":
  parser = OptionParser()
  parser.disable_interspersed_args()
  parser.usage = "usage: %prog [options] command ... "
  parser.add_option("-t", "--test", dest="istest", action="store_true",
      help="show commands to execute ")
  parser.add_option("-d", "--dir", dest="dir", default=".",
      help="root directory to find files recursively")
  parser.add_option("-p", "--pattern", dest="pattern",
      help="file pattern to search")
  (opts, args) = parser.parse_args()
  if (len(args) == 0):
    parser.print_help()
  else:
    cwd = os.getcwd()
    launch(opts, args)
    os.chdir(cwd)
