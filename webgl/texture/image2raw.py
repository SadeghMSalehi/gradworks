#!/usr/bin/python

import Image
from optparse import OptionParser
import io

def main(opt, args):
	im = Image.open(args[0])
	print im.format, im.size, im.mode
	f = io.open(args[1], "wb", 1)
	d = im.getdata()
	for o in d:
		f.write(o)
	f.close()

if (__name__ == "__main__"):
	parser = OptionParser(usage="%prog [options] input output")
	(opt, args) = parser.parse_args()

	if (len(args) == 0):
		parser.print_help()
	else:
		main(opt, args)
