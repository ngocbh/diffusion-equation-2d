"""Solve the 1D diffusion equation using CN and finite differences."""
from time import sleep

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

# The total number of nodes
numx = 5
numy = 5
nnodes = numx*numy
# The total number of times
ntimes = 100
# The time step
dt = 0.5
# The diffusion constant
D = 0.1
# The spatial mesh size
h = 1.0

G = nx.grid_graph(dim=[numx,numy])

# Setup the node list
nodelist = [(i,j) for i in range(numx) for j in range(numy)]

def node_to_pos2d(node):
    """Give a node number, return the point on the physical mesh."""
    return nodelist[node]

def pos_to_node2d(pos):
    """Give the position on the physical mesh, return the node."""
    return nodelist.index(pos)


L = np.matrix(nx.laplacian(G,nodelist=nodelist))

# The rhs of the diffusion equation
rhs = -D*L/h**2

# The main temperature array.
T = np.matrix(np.zeros((nnodes,ntimes)))

# Setting initial temperature. Here we create a temporary array that has the
# physical shape. We then set the temp using x,y coordinates and then reshape
# for the computation. This rehaping matches the nodelist above.
Tnot = np.matrix(np.zeros((numx,numy)))
Tnot[2,2] = 100.0
Tnot.shape = (nnodes,1)
T[:,0] = Tnot

# Set the diffusion constant
D = np.matrix(np.zeros((numx,numy)))
D[:,:] = 1.0 # 1.0 everywhere
D[2,2] = 2.0 # 2.0 somewhere
D.shape = (nnodes,1)


# Setup the time propagator. In this case the rhs is time-independent so we
# can do this once.
ident = np.matrix(np.eye(nnodes,nnodes))
pmat = ident+(dt/2.0)*rhs
mmat = ident-(dt/2.0)*rhs
propagator = np.linalg.inv(mmat)*pmat

# Propagate
for i in range(ntimes-1):
    T[:,i+1] = propagator*T[:,i]
    print "Energy: ", T[:,i].sum()


# To use this to plot call it like this: plot_temp(T,10)

def plot_temp(T, i):
    """Plot the Temperature at the ith time step."""
    Ti = np.copy(T[:,i])
    Ti.shape = (numx,numy)
    plt.contourf(Ti)
    plt.colorbar()
