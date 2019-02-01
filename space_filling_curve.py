import numpy as np
from math import copysign

class FractalCurve:

    def __init__(self,
                 alphabet,
                 proto,
                 base_maps,
                 number_sub,
                 genus = None,
                 dim = None,
                 fractal = None,
                 subdiv_0 = None, 
                 matrix_base_maps = None,
                 subdiv_n = None):
        
        self.alphabet = alphabet
        self.proto = proto
        self.base_maps = base_maps
        self.number_sub = number_sub
        self.genus = self.get_genus()
        self.fractal = self.get_fractal()
        self.dim = self.get_dim()
        self.subdiv_0 = self.get_subdiv_0()
        self.matrix_base_maps = self.get_matrix_base_maps()
        self.subdiv_n = self.get_subdiv_n()

    #Finding the curve genus
    def get_genus(self):
        return len(self.base_maps[0])
    
    #Finding the curve fractality
    def get_fractal(self):
        return len(self.proto)

    #Finding the curve dimension
    def get_dim(self):
        
        #Is the curve mono-fractal?
        index_first = 0 if self.fractal == 1 else 1
        
        #Are fractions reversed in time?
        index_last = 1 if self.base_maps[0][0][-1] in ['0','1'] else 0
            
        return len(self.base_maps[0][0]) - (index_first + index_last)
    
    #Generating the 0-th curve subdivision
    def get_subdiv_0(self):

        #Generating unique vectors matrix
        coordinate_subdiv_0 = np.repeat(np.eye(self.dim),2,axis=0).astype(int)
        for k in range(2*self.dim):
            coordinate_subdiv_0[k,:] = (-1)**k*coordinate_subdiv_0[k,:]
        
        #Generating the 0-th curve subdivision        
        my_dict = {}
        for k in range(len(self.alphabet)):
            my_dict[self.alphabet[k]] = coordinate_subdiv_0[k,:]
    
        subdiv_0 = [np.zeros((self.genus-1,self.dim)) for i in range(self.fractal)] 
        for l in range(self.fractal):
            if len(self.proto[l]) == 1:
                for k in range(self.genus-1):
                    subdiv_0[l][k,:] = my_dict[self.proto[l][0][k]]
            else:
                for k in range(self.genus-1):
                    if len(self.proto[l][k]) == 1:
                        subdiv_0[l][k,:] = my_dict[self.proto[l][k]]
                    else:
                        subdiv_0[l][k,:] = np.sum([my_dict[self.proto[l][k][m]]for m in range(len(self.proto[l][k]))],axis=0)
        return subdiv_0
        
    def get_matrix_base_maps(self):
    
        #Generating coordinates of base maps
        coordinates_base_maps = np.repeat(range(self.dim),2,axis=0).astype(float)        
        for k in range(2*self.dim):
            coordinates_base_maps[k] = (-1)**k*coordinates_base_maps[k]
        
        a0 = 0 if self.fractal == 1 else 1
        
        #Generating matrix base maps
        matrix_base_maps = [np.zeros((self.genus,self.dim+2)) for i in range(len(self.base_maps))]
        for r in range(len(self.base_maps)):
            for k in range(self.genus):
                for m in range(self.dim):
                    for l in range(2*self.dim):
                        if self.base_maps[r][k][m+a0] == self.alphabet[l]:
                            matrix_base_maps[r][k,m+1] = coordinates_base_maps[l]
        
        #Reverse
        for k in range(self.genus):
           if self.fractal == 1:
               if len(self.base_maps[0][0]) == (self.dim+1):
                    matrix_base_maps[0][k,self.dim+1] = float(self.base_maps[0][k][self.dim])
               else:
                    continue
           else:
               for r in range(self.fractal):
                   if len(self.base_maps[0][0]) == (self.dim+2):
                       matrix_base_maps[r][k,self.dim+1] = float(self.base_maps[r][k][self.dim+1])
                   else:
                       continue
        
        #Curve number
        if self.fractal != 1:
            for r in range(self.fractal):
                for k in range(self.genus):
                    matrix_base_maps[r][k,0] = float(self.base_maps[r][k][0])
        return matrix_base_maps

    def get_subdiv_n(self):
        
        subdiv_n = self.subdiv_0
        
        for n in range(self.number_sub):
        
            #Generating fractions
            P = [[[ copysign(1,self.matrix_base_maps[r][k,m+1])*subdiv_n[int(self.matrix_base_maps[r][k,0])][:,int(abs(self.matrix_base_maps[r][k,m+1]))] for k in range(self.genus)] for m in range(self.dim)] for r in range(self.fractal)]
            P = [np.stack(P[r], axis = -1) for r in range(self.fractal)]
        
            #Reversing fractons
            for r in range(self.fractal): 
                for k in range(self.genus):
                    if self.matrix_base_maps[r][k,-1] == 1:
                        P[r][k] = -np.flipud(P[r][k])
                    else:
                        continue
        
            #Glaing fractions
            P1 = [[np.concatenate((P[r][k], [self.subdiv_0[r][k]]), axis = 0) for k in range(self.genus-1)] for r in range(self.fractal)]
            P1 = [np.concatenate(P1[r]) for r in range(self.fractal)]
            P1 = [np.concatenate((P1[r], P[r][-1]), axis = 0) for r in range(self.fractal)]
            
            subdiv_n = P1
            
        #Scaling and shifting the curve
        O = np.zeros((1,self.dim)) 
        subdiv_n = np.concatenate((O, subdiv_n[0]), axis = 0)
        subdiv_n = np.cumsum(subdiv_n,axis = 0)/(self.genus**((self.number_sub+1)/self.dim))
        
        for k in range(self.dim):
            subdiv_n[:,k] = subdiv_n[:,k] - np.amin(subdiv_n[:,k]) + 1/(2*self.genus**((self.number_sub+1)/(self.dim)))
        
        return subdiv_n
