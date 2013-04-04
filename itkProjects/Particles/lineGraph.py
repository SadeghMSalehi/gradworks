#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser


def loadParticles(fIn):
    f = open(fIn)
    l = f.readlines()
    l = [map(float, y.split()[3:5]) for y in l]
    return l


def main(opts, args):
    par = loadParticles(args[0])
    t = np.arange(0, len(par), 1)
    xPos = [x[0] for x in par]
    yPos = [x[1] for x in par]
    plt.plot(t, xPos, 'r--', t, yPos, 'g^', xPos, yPos, 'yo')
    plt.ylabel('position')
    plt.show()
    return


if (__name__ == "__main__"):
    parser = OptionParser("%prog [options] file")
    (opts, args) = parser.parse_args()

    if (len(args) == 0):
        parser.print_help()
    else:
        main(opts, args)
