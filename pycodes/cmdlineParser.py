#!/usr/bin/python

from collections import *
import re
from optparse import OptionParser

pat = re.compile("set[ ]*\(cmdline (.*)\)")

def findOptionName(opt):
  if (opt[0] == '-'):
    if (opt[1] == '-'):
      return opt[2:]
    return opt[1:]
  return opt

def stripVariableName(var):
  if (var[0] == '$'):
    var = var[1:]
  if (var[0] == '{'):
    var = var[1:-1]
  return var

def parseCmd(cmd):
  args = cmd.split(' ')
  executable = args[0]

  print args
  argQ = deque(args[1:])
  
  cmdArgs = []
  while (len(argQ) > 0):
    arg = argQ.popleft()
    if (arg[0] == '-'):
      argOpt = argQ[0]
      if (argOpt[0] != '-'):
        cmdArgs.append((findOptionName(arg), argQ.popleft()))
      else:
        cmdArgs.append((findOptionName(arg), 1))
    else:
      cmdArgs.append((arg,0))
  return (executable, cmdArgs)

def printBMM(appName, cmdArgs):
  f = open("bmm/%s.bmm" % appName, "w")
  str = "<?xml version=1.0?>\n"
  str = str + "<BatchMakeApplicationWrapper>"
  str = str + "\n" + "\t<BatchMakeApplicationWrapperVersion>1.0</BatchMakeApplicationWrapperVersion>"
  str = str + "\n" + "\t<Module>"
  str = str + "\n" + "\t\t<Name>%s</Name>" % appName
  str = str + "\n" + "\t\t<Version>1.0</Version>"
  str = str + "\n" + "\t\t<Path>%s</Path>" % (appName)
  str = str + "\n" + "\t\t<Parameters>"

  argIdx = 0
  for arg in cmdArgs:
    type = 4
    name = arg[0]

    if (arg[1] == 1):
      type = 1
    elif (arg[1] == 0):
      name = "arg%d" % (argIdx)
      argIdx += 1
    else:
      try:
        val = int(arg[1])
        type = 2
      except ValueError:
        try:
          val = float(arg[1])
          type = 3
        except ValueError:
          type = 4

    str = str + "\n" + "\t\t\t<Parameter>"
    if (arg[1] == 0):
      str = str + "\n" + "\t\t\t\t<Type>4</Type>"
      str = str + "\n" + "\t\t\t\t<Name>%s</Name>" % (name)
      str = str + "\n" + "\t\t\t\t<Value>--%s</Value>" % (name)
      str = str + "\n" + "\t\t\t\t<External>%s</External>" % (0) 
      str = str + "\n" + "\t\t\t\t<Optional>%s</Optional>" % (0) 
    else:
      str = str + "\n" + "\t\t\t\t<Type>1</Type>"
      str = str + "\n" + "\t\t\t\t<Name>%s</Name>" % (name)
      str = str + "\n" + "\t\t\t\t<Value>--%s</Value>" % (name)
      str = str + "\n" + "\t\t\t\t<External>%s</External>" % (0) 
      str = str + "\n" + "\t\t\t\t<Optional>%s</Optional>" % (0) 
      if (type > 1):
        str = str + "\n" + "\t\t\t\t<SubParameter>"
        str = str + "\n" + "\t\t\t\t\t<Type>%d</Type>" % (type)
        str = str + "\n" + "\t\t\t\t\t<Name>%svalue</Name>" % (name)
        str = str + "\n" + "\t\t\t\t\t<External>0</External>"
        str = str + "\n" + "\t\t\t\t\t<Optional>0</Optional>"
        str = str + "\n" + "\t\t\t\t</SubParameter>"
    str = str + "\n" + "\t\t\t</Parameter>"
  str = str + "\n" + "\t\t</Parameters>"
  str = str + "\n" + "\t</Module>"
  str = str + "\n" + "</BatchMakeApplicationWrapper>"
  f.write(str)
  f.close()

def printAppOptions(appName, args):
  appName = stripVariableName(appName)
  strOut = ""
  appCanonicalName = "App" + appName[0].upper() + appName[1:]
  strOut = strOut + "\n" + "SetApp(%s @%s)" % (appName, appCanonicalName)
  argIdx = 0
  for arg in args:
    if (arg[1] == 0):
      strOut = strOut + "\n" + "SetAppOption(%s.arg%d %s)" % (appName, argIdx, arg[0])
      argIdx = argIdx + 1
    else:
      if (arg[1] == 1):
        strOut = strOut + "\n" + "SetAppOption(%s.%s 1)" % (appName, arg[0])
      else:
        strOut = strOut + "\n" + "SetAppOption(%s.%s 1)" % (appName, arg[0])
        strOut = strOut + "\n" + "SetAppOption(%s.%s.%s %s)" % (appName, arg[0], arg[0] + "value", arg[1])

  strOut = strOut + "\n" + "Set (cmdAppName %s)" % (appName)
  printBMM(appCanonicalName, args)
  return strOut

def main(opts, argv):
  newCode = ""
  f = open(argv[0])
  for l in f:
    l = l.strip()
    m = pat.match(l)
    if (m):
      newCode = newCode + "\n" + "#" + l
      cmd = m.group(1)
      print "Parsing ... " + cmd
      (exe, opt) = parseCmd(cmd)
      print opt
      appsRun = printAppOptions(exe, opt)
      newCode = newCode + "\n" + appsRun
    else:
      newCode = newCode + "\n" + l
      
  fout = open("out/%s" % argv[0], "w")
  fout.write(newCode)
  fout.close()

  return

parser = OptionParser(usage="%prog [options] file1 file2 ...")
(opts, argv) = parser.parse_args()

main([], argv)
