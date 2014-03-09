"NIRAL lab python utility package"
#!/usr/bin/python

import os, os.path, subprocess, sys
import fnmatch 

def colorize(text, color_code): return color_code + text + "\033[0m" 
def red(t): return colorize(t, "\033[1m\033[31m")
def green(t): return colorize(t, "\033[1m\033[32m")
def darkgreen(t): return colorize(t, "\033[32m")
def yellow(t): return colorize(t, "\033[1m\033[33m")
def blue(t): return colorize(t, "\033[1m\033[34m")
def darkblue(t): return colorize(t, "\033[34m")
def pur(t): return colorize(t, "\033[1m\033[35m")

def findFilesByName(dir, fargs):
  "find all files under the dir with given name list"
  freturns = []
  for root, dir, files in os.walk(dir):
    for f in fargs:
      if f in files:
        freturns.append(os.path.join(root,f))

  freturns.sort()
  return freturns

def findFilesByPattern(dir, fpattern):
  "find all files with the pattern from dir"
  freturns = []
  for root, dir, files in os.walk(dir):
    freturns.extend(map(lambda x: os.path.join(root,x), fnmatch.filter(files, fpattern)))

  freturns.sort()
  return freturns

def exe(cmd, args):
  try:
    cmdline = cmd + " " + args
    print yellow(cmd) + " " + args

    retcode = subprocess.call(cmdline, shell=True)
    if retcode < 0:
      print >>sys.stderr, "Child was terminated by signal", -retcode
      return 1
    else:
      return 0
  except OSError, e:
        print >>sys.stderr, "Execution failed:", e
        return 2
