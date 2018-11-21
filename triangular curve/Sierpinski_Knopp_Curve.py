import sys
import os
sys.path.append(os.path.dirname(sys.argv[0])+'/..')

import numpy as np
from pl_curves import plcurves

#Indicate subdivision number and run
k = 4

#The definition of basis vectors and the origin
O = np.array([[0,0]])
i = np.array([[1,0]])
j = np.array([[0,1]])
ij = np.array([[1,1]])
iJ = np.array([[1,-1]])

#Curve Prototype (null curve subdivision)
subdiv_0 = np.concatenate([j, ij, i])

#Generating the n-th subdivision of the curve
subdiv_n = subdiv_0
for n in range(k):
    
    #Base mapping of second fraction
    Ji = np.column_stack((-subdiv_n[:,1],subdiv_n[:,0]))

    #Base mapping of third fraction
    jI = np.column_stack((subdiv_n[:,1],-subdiv_n[:,0]))
    
    #Generating following subdivision
    subdiv_n = np.concatenate([subdiv_n, j, Ji, ij, jI, i, subdiv_n])

#Joining the two triangles
subdiv_n = np.concatenate([subdiv_n,iJ,-subdiv_n])

#Generating, scaling and shifting the curve coordinates
subdiv_n = np.cumsum(np.concatenate([O,subdiv_n]),axis = 0)/2**(k+2)
subdiv_n = np.column_stack((subdiv_n[:,0],subdiv_n[:,1] + subdiv_n[1,1])) + 1/2**(k+3)

plcurves(subdiv_n, 2, 4, k)