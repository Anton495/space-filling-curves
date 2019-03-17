class FractalCurve:
    
    def __init__(self,chain_proto,base_maps,
                 alphabet=None,dim=None,genus=None,fractal=None):
        
        self.chain_proto = chain_proto
        self.base_maps = base_maps
        self.alphabet = alphabet if alphabet is not None else self.get_alphabet()
        self.dim = dim if dim is not None else self.get_dim()
        self.genus = genus if genus is not None else self.get_genus()
        self.fractal = fractal if fractal is not None else self.get_fractal()

    def get_alphabet(self):
        return 'ijklmnop'

    def get_dim(self):
        '''get the curve dimension'''
        return len(set(''.join(self.chain_proto[0]).lower()))

    def get_genus(self):
        '''get the curve genus'''
        return len(self.base_maps[0])

    def get_fractal(self):
        '''get the curve fractality'''
        return len(self.chain_proto)

    def get_curve_coord(self,chain_code):
        '''chain code => vectors => coordinates'''
        
        # Формируем словарь единичных векторов
        vect_dict = {}
        for k in range(self.dim):
            coord = [0] * self.dim
            coord[k] = 1
            vect_dict[self.alphabet[k]] = coord
            vect_dict[self.alphabet[k].upper()] = [-m for m in coord]
        
        # Эта функция получает диагональный вектор из нескольких единичных векторов
        # путем суммирования их координат. Например ij = i + j = [1,0] + [0,1] = [1,1]
        def get_diag_coord(vector):
            arg = [vect_dict[k] for k in vector]
            coord = list(map(sum,zip(*arg)))
            return coord
        
        # Переходим от цепного кода к единичным векторам (применяем словарь)
        a = len(''.join([''.join(k) for k in self.chain_proto]))
        if a//self.fractal == a/self.fractal:
            subdiv = list(map(vect_dict.get,chain_code))
        else:
            subdiv = [get_diag_coord(m) for m in chain_code]
        
        # Переходим от единичных векторов к координатам кривой (суммируем координаты)
        curve_coord = [[0]*self.dim] + subdiv        
        for l in range(len(curve_coord)-1):
            curve_coord[l+1] = [c + d for c, d in zip(curve_coord[l], curve_coord[l+1])]
        
        return curve_coord
    
    def get_proto(self):
        '''get the curve prototype'''
        return [self.get_curve_coord(k) for k in self.chain_proto]
    
    def get_div(self):
        '''get the curve div'''
        return max(self.get_curve_coord(self.chain_proto[0]))[0]+1
    
    def get_dict_bm(self,id_bm,bm):
        '''get base map dictonary'''    
        dict_bm={}
        for k in range(self.dim):
            m = bm.lower().index(id_bm[k])
            letter = id_bm[m]
            letter = letter if id_bm[k] == bm[m] else letter.upper()
            dict_bm[id_bm[k]] = letter
            dict_bm[id_bm[k].upper()] = letter.swapcase()
                        
        return dict_bm
    
    def get_fraction(self,sub,dic,inv):
        '''apply base map and reverse to some curve fraction'''
        fraction = [''.join(list(map(dic.get, k))) for k in sub]
        
        if inv == '1':
            fraction = list(reversed(fraction))
            fraction = [k.swapcase() for k in fraction]
        
        return fraction
    
    def get_subdiv(self,sub_numb):
        '''get n-th curve subdivision'''
        
        # Добавляем номер кривой в базовые преобразования для монофракталов
        if self.fractal != 1:   
            base_maps = self.base_maps
        else:
            base_maps = [['0' + k for k in self.base_maps[0]]]  
        
        # Формируем список словарей для всех базовых преобразований
        id_bm = self.alphabet[:self.dim]
        list_dict = [[self.get_dict_bm(id_bm,base_maps[k][m][1:]) 
                      for m in range(self.genus)] for k in range(self.fractal)]        
        
        # Определяем нулевое подразделение кривой
        sub_k = self.chain_proto
        
        for n in range(sub_numb):
            
            # Формируем список преобразованных фракций
            sub_n = [[self.get_fraction(sub_k[int(base_maps[k][m][0])],list_dict[k][m],base_maps[k][m][-1]) 
                      for m in range(self.genus)] for k in range(self.fractal)]
            
            # Добавляем связующие ребра между фракциями
            [[sub_n[k].insert(2*m+1,[self.chain_proto[k][m]])
              for m in range(self.genus-1)] for k in range(self.fractal)]
            
            # Объединяем фракции и связующие ребра в один список
            sub_n = [sum(k,[]) for k in sub_n]
            
            # Переопределяем (n-1)-ое на n-ое подразделение
            sub_k = sub_n.copy()
            
        subdiv = self.get_curve_coord(sub_k[0])
        
        return subdiv
