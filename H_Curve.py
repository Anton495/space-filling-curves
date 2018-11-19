import numpy as np
from pl_curves import plcurves

#Indicate subdivision number and run
k = 4

#The definition of basis vectors and the origin
O = np.array([[ 0 , 0 ]])
i = np.array([[ 1 , 0 ]])
j = np.array([[ 0 , 1 ]])

#Curve Prototype (null curve subdivision)
subdiv_0 = j

#Generating the n-th subdivision of the curve
subdiv_n = subdiv_0
for n in range(k):
    
    #Base mapping of second fraction and reverse
    Ij = np.column_stack((-subdiv_n[:,0],subdiv_n[:,1]))
    Ij1 = np.flipud(Ij)

    #Base mapping of third fraction and reverse
    iJ = np.column_stack((subdiv_n[:,0],-subdiv_n[:,1]))
    iJ1 = np.flipud(iJ)
    
    #Generating following subdivision
    subdiv_n = np.concatenate([subdiv_n, j, Ij1, i, iJ1, i, subdiv_n])

#Joining the two triangles
subdiv_n = np.concatenate([subdiv_n,i,-subdiv_n])

#Generating, scaling and shifting the curve coordinates
subdiv_n = np.cumsum(np.concatenate([O,subdiv_n]),axis = 0)/2**(k+1) + 1/2**(k+2)

plcurves(subdiv_n, 2, 4, k)