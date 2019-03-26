from numpy import cumsum
import itertools as it

class FractalCurve:
    
    def __init__(self,chain_proto,base_maps,
                 alphabet=None,dim=None,genus=None,fractal=None,vect_dict=None):
        
        self.chain_proto = chain_proto
        self.base_maps = base_maps
        self.alphabet = alphabet if alphabet is not None else self.get_alphabet()
        self.dim = dim if dim is not None else self.get_dim()
        self.genus = genus if genus is not None else self.get_genus()
        self.fractal = fractal if fractal is not None else self.get_fractal()
        self.vect_dict = vect_dict if vect_dict is not None else self.get_vect_dict()

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
        fractal = len(self.chain_proto)
        # Добавляем номер кривой в базовые преобразования для монофракталов
        if fractal == 1:
            self.base_maps = [['0' + k for k in self.base_maps[0]]]  
        return fractal

    def get_vect_dict(self):
        '''get unit and diagonal vectors dictonary'''
        
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
        
        # Проверяем наличие диагональных шагов в chain_proto
        a = len(''.join([''.join(k) for k in self.chain_proto]))
        if a//self.fractal != a/self.fractal:
            # Добавляем диагональные шаги в словарь
            all_letters = [[self.alphabet[k],self.alphabet[k].upper()] for k in range(self.dim)]
            for m in range(2,self.dim+1):
                all_comb_let = list(it.permutations(all_letters,m))
                all_diag_coord = list(map(''.join, it.chain(*[it.product(*k) for k in all_comb_let])))
                for k in all_diag_coord:
                    vect_dict[k] = get_diag_coord(k)
         
        return vect_dict
            
    def get_curve_coord(self,chain_code,start=None):
        '''chain code => vectors => coordinates'''
        
        # Переходим от цепного кода к единичным векторам (применяем словарь)    
        subdiv = list(map(self.vect_dict.get,chain_code))
        
        # Определяем начальную координату для кривой
        if start == None:
            curve_coord = [[0]*self.dim] + subdiv
        else:
            curve_coord = [start] + subdiv
            
        # Переходим от единичных векторов к координатам кривой (суммируем координаты)
        # 1 способ
        #curve_coord = list(zip(*map(it.accumulate, zip(*curve_coord))))
        # 2 способ
        # Коммутативная сумма быстрее на 20%-30%, но возвращает не list of tuple, а numpy.array
        curve_coord = cumsum(curve_coord,axis = 0)
        
        return curve_coord
    
    def get_proto(self):
        '''get the curve prototype'''
        return [self.get_curve_coord(k) for k in self.chain_proto]
    
    def get_div(self):
        '''get the curve div'''
        return max(self.get_curve_coord(self.chain_proto[0]))[0]+1
    
    def get_fraction(self,sub,bm):
        '''apply base map and reverse to some curve fraction'''
        
        # Определяем тождественое базовое преобразование
        id_bm = self.alphabet[:self.dim]
        
        # Создаем словарь базового преобразования
        dict_bm={}
        for k in range(self.dim):
            m = bm.lower().index(id_bm[k])
            letter = id_bm[m]
            letter = letter if id_bm[k] == bm[m] else letter.upper()
            dict_bm[id_bm[k]] = letter
            dict_bm[id_bm[k].upper()] = letter.swapcase()
            
        # Поворачиваем фракцию
        fraction = [''.join(list(map(dict_bm.get, k))) for k in sub]
        
        # Обращаем по времени
        if bm[-1] == '1':
            fraction = list(reversed(fraction))
            fraction = [k.swapcase() for k in fraction]
        
        return fraction
    
    def get_subdiv(self,sub_numb,plot=True):
        '''get n-th curve subdivision'''
        
        # Задаем нулевое подразделение кривой
        sub_k = self.chain_proto
        for n in range(sub_numb):

            # Формируем список преобразованных фракций
            sub_n = [[self.get_fraction(sub_k[int(self.base_maps[k][m][0])],self.base_maps[k][m][1:]) 
                      for m in range(self.genus)] for k in range(self.fractal)]
            
            # Добавляем связующие ребра между фракциями для графика
            if plot==True:    
                [[sub_n[k].insert(2*m+1,[self.chain_proto[k][m]])
                  for m in range(self.genus-1)] for k in range(self.fractal)]

            # Объединяем все фракции в одну фракцию
            sub_n = [sum(k,[]) for k in sub_n]
            
            # Определяем (n-1)-ое как n-ое подразделение
            sub_k = sub_n.copy()
            
        subdiv = self.get_curve_coord(sub_k[0])
        
        return subdiv
