import MDAnalysis as mda
import numpy as np
from sys import argv

if not (argv[1] or argv[2]):
    print "usage: python nowater.py test.xplor.psf test.dcd"
    exit()
else:
    psf = str(argv[1])
    dcd = str(argv[2])
    jobname = dcd[:-4]

u = mda.Universe(psf,dcd)

nowater = u.select_atoms('all and not resname SWM4')

with mda.Writer("%s.nowater.dcd" % (jobname), nowater.n_atoms) as W:
    for ts in u.trajectory:
        W.write(nowater)
