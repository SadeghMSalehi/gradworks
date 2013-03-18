#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
import sys, os


def bargraphAB(opts,args):
	l = sys.stdin.readlines()
	l = [ r.strip().split(" ") for r in l ]
	d = [map(lambda x: float(x), s[1:]) for s in l]
	n = [ x[0] for x in l]
	t = [ int(x[3]) for x in d]

	a = [ 100*x[0]/z for (x,z) in zip(d,t)]
	ab = [ 100*x[1]/z for (x,z) in zip(d,t)]
	b = [ 100*x[2]/z for (x,z) in zip(d,t)]

	a_bottom = [ x+y for (x,y) in zip(ab,b) ]
	ab_bottom = b
	ind = np.arange(len(l))
	width = 0.35
	p1 = plt.bar(ind, a, width, color='c', bottom=a_bottom)
	p2 = plt.bar(ind, ab, width, color='y', bottom=ab_bottom)
	p3 = plt.bar(ind, b, width, color='m')

	mm = 0
	mx = 100

	plt.ylabel('Similarity Scores')
	plt.title('Scores of FP/TP/FN')
	plt.xticks(ind+width/2., n)
	plt.yticks(np.arange(mm,mx,(mx-mm)/10))
	plt.legend( (p1[0], p2[0], p3[0]), ('FP', 'TP', 'FN') )

	plt.savefig("foo.pdf")

def main(opts, args):
	bargraphAB(opts,args)

if (__name__ == "__main__"):
	parser = OptionParser(usage="%prog [options] ")
	(opts, args) = parser.parse_args()
	main(opts, args)
"""
N = 5
menMeans   = (20, 35, 30, 35, 27)
womenMeans = (25, 32, 34, 20, 25)
menStd     = (2, 3, 4, 1, 2)
womenStd   = (3, 5, 2, 3, 3)
ind = np.arange(N)    # the x locations for the groups
width = 0.35       # the width of the bars: can also be len(x) sequence

p1 = plt.bar(ind, menMeans,   width, color='r', yerr=womenStd)
p2 = plt.bar(ind, womenMeans, width, color='y',
             bottom=menMeans, yerr=menStd)

plt.ylabel('Scores')
plt.title('Scores by group and gender')
plt.xticks(ind+width/2., ('G1', 'G2', 'G3', 'G4', 'G5') )
plt.yticks(np.arange(0,81,10))
plt.legend( (p1[0], p2[0]), ('Men', 'Women') )

plt.show()
"""
