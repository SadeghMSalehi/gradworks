#!/usr/bin/python

import os,sys,glob
from optparse import OptionParser

def main(opts,args):
	files = args
	fileHandles = []
	for f in files:
		fo = open(f,"r")
		fileHandles.append(fo)
	out = open(opts.output, "w")
	isEnd = False
	while not isEnd:
		buf = ""
		isEnd = True
		for fh in fileHandles:
			line = fh.readline()
			line = line.strip()
			if (line == "" or line.startswith("#")):
				data = opts.emptyColumn
			else:
				cols = line.split(opts.delimiter)
				if (len(cols) > opts.colIndex):
					data = cols[opts.colIndex]
					data = data.strip()
				buf = buf + data + "\t"
				isEnd = False
		out.write(buf)
		out.write("\n")
	out.close()

if (__name__ == "__main__"):
	parser = OptionParser(usage="""%prog [options] 1.txt 2.txt ...
		Extract a column specified by -c or --columnIndex and append vertically in the order of args 	
		For example, 
			> cat 1.txt
				1  2  3
				2  4  1
				3  5
			> cat 2.txt
				3  5  7
				4  9  3
	    > %prog -o out.txt -c 1 -e "" 1.txt 2.txt 
			> cat out.txt
				2 5
				4 9
				5 
		If there is no column exists, or one file is longer than another, the empty column data will 
		filled as the value of -e or --emptyColumn
	""")
	parser.add_option("-o", "--output", help="specify an output file name", dest="output", default="output.txt")
	parser.add_option("-c", "--columnIndex", help="specify an index of the column to extract", dest="colIndex", default=1, type="int")
	parser.add_option("-e", "--emptyColumn", help="specify a string to fill in empty column", dest="emptyColumn", default="")
	parser.add_option("-d", "--delimiter", help="specify a delimiter to split into columns", dest="delimiter", default=" ")
	(opts,args) = parser.parse_args()
	if (len(args) == 0):
		parser.print_help()
	else:
		main(opts,args)
