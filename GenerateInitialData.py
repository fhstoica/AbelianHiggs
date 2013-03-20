#! /usr/bin/python

grid_points = 64

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
            f.write("%s %s %s %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n" % (a, b, c,
                                                                        (random.random()-0.5),
                                                                        (random.random()-0.5),
                                                                        (random.random()-0.5),
                                                                        (random.random()-0.5),
                                                                        (random.random()-0.5),
                                                                        (random.random()-0.5))
                    )
f.close()
