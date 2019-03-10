import numpy as np
from math import copysign

class FractalCurve:

    def __init__(self,chain_code,base_maps,
                 genus = None,dim = None,fractal = None,base_maps_matrix = None):

        self.chain_code = chain_code
        self.base_maps = base_maps
        self.dim = dim if dim is not None else self.get_dim()
        self.genus = genus if genus is not None else self.get_genus()
        self.fractal = fractal if fractal is not None else self.get_fractal()
        self.base_maps_matrix = base_maps_matrix if base_maps_matrix is not None else self.get_base_maps_matrix()

    #Finding the curve genus
    def get_genus(self):
        return len(self.base_maps[0])
    
    #Finding the curve fractality
    def get_fractal(self):
        return len(self.chain_code)

    #Finding the curve dimension
    def get_dim(self):
        return len(set(''.join(self.chain_code[0]).lower()))
    
    #Generating the 0-th curve subdivision
    def get_subdiv_0(self):

        alphabet = 'ijklmn'
        assert self.dim <= 6
        
        A = [0]*self.dim
        A = [list(A) for k in A]
        for k in range(self.dim):
            A[k][k] = 1

        B = [0]*self.dim
        B = [list(B) for k in B]
        for k in range(self.dim):
            B[k][k] = -1

        my_dict = {}
        for k in range(self.dim):    
            my_dict[alphabet[k]] = A[k]

        for k in range(self.dim):
            my_dict[alphabet.upper()[k]] = B[k]
        
        def diag_coord(vector):
            C = [my_dict[k] for k in vector]
            coord = list(map(sum,zip(*C)))
            return coord
        
        subdiv_0 = [list(map(my_dict.get, self.chain_code[k][0])) if len(self.chain_code[k]) == 1 
                 else [diag_coord(m) for m in self.chain_code[k]] for k in range(self.fractal)]
        
        return subdiv_0
    
    #Generating the base maps matrix
    def get_base_maps_matrix(self):
    
        alphabet = 'ijklmn'
        assert self.dim <= 6
        
        my_dict = {}
        for k in range(self.dim):    
            my_dict[alphabet[k]] = float(k)
        
        for k in range(self.dim):
            my_dict[alphabet.upper()[k]] = -float(k)
            
        for k in range(self.fractal+1):
            my_dict[str(k)] = k
            
        base_maps_matrix = [[list(map(my_dict.get, self.base_maps[m][k])) 
                            for k in range(self.genus)] for m in range(self.fractal)]
        
        [k.insert(0,0) for k in base_maps_matrix[0] if self.fractal == 1]
        
        for r in range(self.fractal):
            for k in range(self.genus):
                if self.base_maps[r][k][-1] not in ['0','1']:
                    base_maps_matrix[r][k].append(0)
                    
        return base_maps_matrix
    
    #Generating the n-th curve subdivision
    def get_subdiv_n(self,subdiv_number):
        
        subdiv_n = np.array(self.get_subdiv_0())
        
        for n in range(subdiv_number):
        
            #Generating fractions
            P = [[[ copysign(1,self.base_maps_matrix[r][k][m+1])*subdiv_n[int(self.base_maps_matrix[r][k][0])][:,int(abs(self.base_maps_matrix[r][k][m+1]))] 
                for k in range(self.genus)] for m in range(self.dim)] for r in range(self.fractal)]
            P = [np.stack(P[r], axis = -1) for r in range(self.fractal)]
        
            #Reversing fractons
            for r in range(self.fractal): 
                for k in range(self.genus):
                    if self.base_maps_matrix[r][k][-1] == 1:
                        P[r][k] = -np.flipud(P[r][k])
                    else:
                        continue
        
            #Glaing fractions
            P1 = [[np.concatenate((P[r][k], [self.get_subdiv_0()[r][k]]), axis = 0) for k in range(self.genus-1)] for r in range(self.fractal)]
            P1 = [np.concatenate(P1[r]) for r in range(self.fractal)]
            P1 = [np.concatenate((P1[r], P[r][-1]), axis = 0) for r in range(self.fractal)]
            
            subdiv_n = P1
        
        #Scaling and shifting the curve
        O = np.zeros((1,self.dim)) 
        subdiv_n = np.concatenate((O, subdiv_n[0]), axis = 0)
        subdiv_n = np.cumsum(subdiv_n,axis = 0)/(self.genus**((subdiv_number+1)/self.dim))
        
        for k in range(self.dim):
            subdiv_n[:,k] = subdiv_n[:,k] - np.amin(subdiv_n[:,k]) + 1/(2*self.genus**((subdiv_number+1)/(self.dim)))
        return subdiv_n
