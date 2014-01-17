#!/usr/bin/python

import os,sys

if (__name__ == "__main__"):
	if (len(sys.argv) < 2):
		sys.exit(0)
	n = sys.argv[1]
	o = sys.argv[2]
	f = open(n)	
	g = open(o, "w")
	l = f.readlines()
	isMD = False
	for r in l:
		if (r.startswith("/*MD")):
			isMD = True
		elif (r.startswith("*/") and isMD):
			isMD = False
		elif (isMD):
			g.write(r[1:])
	f.close()
	g.close()

