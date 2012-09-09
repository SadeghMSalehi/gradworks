#!/usr/bin/python

import fileinput, glob, string, sys, os
from os.path import join
from optparse import OptionParser
import niralutil as niral

# replace a string in multiple files

def main(opts, args):
  stext = args[0]
  rtext = args[1]
  pattern = opts.pattern

  print "finding [" + pattern + "]: '" + stext + "' replacing with: '" + rtext + "' in: "

  if (opts.r):
    print "recursive search..."
    files = niral.findFilesByPattern(".", pattern)
  else:
    files = glob.glob(pattern)

  if (len(files) == 0):
    print "do nothing..."
    return

  print "\n".join(files)

  for line in fileinput.input(files, inplace = 1):
    lineno = 0
    lineno = string.find(line, stext)
    if lineno > 0:
      line = line.replace(stext, rtext)
    sys.stdout.write(line)

if (__name__ == "__main__"):
  parser = OptionParser(usage="%prog [options] search_term replace_term")
  parser.add_option("-r", "--recursive", dest="r", help="enable recursive search of the pattern", action="store_true", default=False)
  parser.add_option("-p", "--pattern", dest="pattern", help="file pattern in which strings are replaced", default="*")
  (opts, args) = parser.parse_args()
  if (len(args) < 2):
    parser.print_help()
  else:
    main(opts, args)
