#!/usr/bin/python3


import numpy as np

import matplotlib
# matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib import rc
rc("text", usetex=True)
rc("font",**{"family":"serif","serif":["Times"]})
rc("font",**{"family":"sans-serif","sans-serif":["Helvetica"]})

import sys


e1 = [int(i) for i in np.loadtxt("meshfiles/edg1.dat")]
e2 = [int(i) for i in np.loadtxt("meshfiles/edg2.dat")]

fig, ax = plt.subplots()

for i in [20]: # range(1,26):
    
    print(i)

    xy = np.loadtxt("outpt/out" + str(i) + ".m")
    x, y = xy[0,], xy[1,]

    ax.clear()
    for j in range(len(e1)):
        ax.plot([x[e1[j]-1],x[e2[j]-1]], [y[e1[j]-1],y[e2[j]-1]], color="b", linewidth=0.5)

    ax.set_aspect("equal")
    plt.savefig("outpt/out" + str(i) + ".pdf", bbox_inches="tight")
