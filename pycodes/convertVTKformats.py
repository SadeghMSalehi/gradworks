#!/biomed-resimg/crews_rodent/devel/linux/ParaView3/ParaView-3.10.0-Linux-x86_64/bin/pvpython
import niral
import sys

p = niral.readVTK(sys.argv[1])
niral.writeVTK(sys.argv[2], p)
