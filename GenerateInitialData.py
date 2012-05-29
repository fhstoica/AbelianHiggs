#! /usr/bin/python

grid_points     = 64

import sys
import os
import random

OUTFILE      = 'initial_data.dat'
r_seed       = float(sys.argv[1])

random.seed(r_seed)
f  = open(OUTFILE , 'w')
for a in range(grid_points) :
    for b in range(grid_points) :
        for c in range(grid_points) :
            f.write(str(a)+' '+str(b)+' '+str(c)+' ' 
                    +str("%1.4f" % (random.random() - 0.5))+' '
                    +str("%1.4f" % (random.random() - 0.5))+'\n')
            continue
        continue
    continue
f.close()
