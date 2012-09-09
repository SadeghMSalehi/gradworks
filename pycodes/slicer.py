#!/usr/bin/python

import os,sys
import niralutil as nu

slicerPath="/tools/Slicer3/Slicer3-3.6-2010-06-10-linux-x86_64/"

args = " ".join(sys.argv[1:])
nu.exe(slicerPath + "Slicer3 ", args)
