#!/usr/bin/python

import sys

def readLines(fname):
  f = open(fname, "r")
  lines = map(lambda x: x.rstrip(), f.readlines())
  return lines

if (__name__ == "__main__"):
  contents = []

  for i in range(1, len(sys.argv)):
    lines = readLines(sys.argv[i]);
    if (i > 1):
      if (len(contents[0]) != len(lines)):
        print "Different number of lines (" + sys.argv[i] + " ): " + str(len(contents[0])) + " != " + str(len(lines))
        sys.exit(-1)
    contents.append(readLines(sys.argv[i]))

  for i in range(0, len(contents[0])):
    for j in range(0, len(contents)):
      sys.stdout.write(contents[j][i])
      sys.stdout.write("\t")
    sys.stdout.write("\n")
