import itertools as it
from functools import reduce
from math import gcd
import sympy as sp

class FractalCurve:
    
    def __init__(self,chain_proto,base_maps,
                 dim=None,alph=None,genus=None,div=None,fractal=None,vect_dict=None):
        
        self.chain_proto = chain_proto
        self.base_maps = base_maps
        self.dim = dim if dim is not None else self.get_dim()
        self.alph = alph if alph is not None else self.get_alphabet()
        self.genus = genus if genus is not None else self.get_genus()
        self.div = div if div is not None else self.get_div()
        self.fractal = fractal if fractal is not None else self.get_fractal()
        self.vect_dict = vect_dict if vect_dict is not None else self.get_vect_dict()

    def get_dim(self):
        '''get the curve dimension'''
        return len(set(''.join(self.chain_proto[0]).lower()))

    def get_alphabet(self):
        '''get alphabet for curve prototypes and base maps'''
        return 'ijklmnop'[:self.dim]

    def get_genus(self):
        '''get the curve genus'''
        return len(self.base_maps[0])

    def get_div(self):
        '''get the curve div'''
        return int(self.genus**(1/self.dim)) #max(self.get_curve_coord(self.chain_proto[0]))[0]+1

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
            vect_dict[self.alph[k]] = coord
            vect_dict[self.alph[k].upper()] = [-m for m in coord]
        
        # Эта функция получает диагональный вектор из нескольких единичных векторов
        # путем суммирования их координат. Например ij = i + j = [1,0] + [0,1] = [1,1]
        def get_diag_coord(vector):
            arg = [vect_dict[k] for k in vector]
            coord = list(map(sum,zip(*arg)))
            return coord
        
        # Добавляем диагональные шаги в словарь
        all_letters = [[self.alph[k],self.alph[k].upper()] for k in range(self.dim)]
        for m in range(2,self.dim+1):
            all_comb_let = list(it.permutations(all_letters,m))
            all_diag_coord = list(map(''.join, it.chain(*[it.product(*k) for k in all_comb_let])))
            for k in all_diag_coord:
                vect_dict[k] = get_diag_coord(k)
        
        return vect_dict

    def get_proto(self):
        '''get the curve prototypes'''
        return [self.get_curve_coord(k) for k in self.chain_proto]
    
    def get_curve_coord(self,chain_code,start=None):
        '''chain code => vectors => coordinates'''
        # Переходим от цепного кода к единичным векторам (применяем словарь)    
        subdiv = list(map(self.vect_dict.get,chain_code))
                    
        # Определяем начальную координату для кривой
        curve_coord = [[0]*self.dim] + subdiv if start == None else [start] + subdiv
        
        # Переходим от единичных векторов к координатам кривой (суммируем координаты)
        curve_coord = list(zip(*map(it.accumulate, zip(*curve_coord))))
        
        return curve_coord
    
    def get_fraction(self,sub,bm):
        '''apply base map and reverse to some curve fraction'''
        # Создаем словарь базового преобразования, например kIJ = {i:J,j:K,k:i}
        dict_bm={}
        for k in range(self.dim):
            dict_bm[bm[k]] = self.alph[k]
            dict_bm[bm[k].swapcase()] = self.alph[k].upper()
            
        # Поворачиваем фракцию (применяем словарь)
        fraction = [''.join(list(map(dict_bm.get, k))) for k in sub]
        
        # Обращаем по времени
        if bm[-1] == '1':
            fraction = list(reversed(fraction))
            fraction = [k.swapcase() for k in fraction]
        
        return fraction
    
    def vertex_chain_proto(self):
        '''Функция находит цепной код, построенный по вершинам прототипов'''
        # Находим вершины в прототипах
        all_vertices = self.get_all_vertices()
        
        # Переходим от координат к векторам 
        all_vectors = []
        for r in range(self.fractal):
            vectors = []
            for l in range(2**self.dim-1):
                vectors.append([x - y for x, y in zip(all_vertices[r][l+1], all_vertices[r][l])])
            all_vectors.append(vectors)
        
        # Масшабируем вектора (приводим в соответствие со словарем)
        scaling_vector = [self.div-1]*self.dim
        for r in range(self.fractal): 
            for l in range(2**self.dim-1):    
                all_vectors[r][l] = str([x//y for x, y in zip(all_vectors[r][l], scaling_vector)])

        # Формируем обратный словарь
        inv_vect_dict = {str(v): k for k, v in self.vect_dict.items()}

        return [list(map(inv_vect_dict.get,r)) for r in all_vectors]
    
    def get_subdiv(self,sub_numb,plot=True):
        '''get n-th curve subdivision'''
        # Определяем нулевое подразделение кривой
        sub_k = self.chain_proto if plot==True else self.vertex_chain_proto()
        
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
            
            # Определяем (n-1)-ое подразделение как n-ое подразделение
            sub_k = sub_n.copy()
        
        return sub_k
    
    def get_con_junction(self,id_bm,junction):
        '''Функция приводит стык к каноническому виду'''
        bm2 = ''
        for l in range(1,self.dim+1):
            m = junction[0].lower().index(id_bm[l])
            if junction[0][m] == junction[0][m].upper():
                bm2 = bm2 + junction[1][m].swapcase()
            else:
                bm2 = bm2 + junction[1][m]
                
        if junction[0][-1] == '1': id_bm = id_bm + '1'
        if junction[1][-1] == '1': bm2 = bm2 + '1'
            
        bm1 = junction[0][0] + id_bm[1:]
        bm2 = junction[1][0] + bm2
            
        con_junction = (bm1,bm2)
        
        return con_junction

    def get_bm_bm(self,id_bm,bm,bm2):
        '''Возвращает базовое преобразование базового преобразования. Например bm2(bm) = KiJ(jkI)'''
        new_bm = ''
        for l in range(1,self.dim+1):
            m = id_bm.index(bm2[l].lower())   
            if bm[m]  == bm[m] .upper():
                new_bm = new_bm + bm[m].upper() if bm2[l] == bm2[l].lower() else new_bm + bm[m].lower()
            else:
                new_bm = new_bm + bm[m].upper() if bm2[l] == bm2[l].upper() else new_bm + bm[m].lower() 
        
        if bm[-1] == '1': new_bm = new_bm + '1'
            
        new_bm = bm[0] + new_bm
        
        return new_bm

    def get_bms(self,base):
        '''Функция находит базовые преобразования производных стыков'''
        new_base = []
        for r in range(self.fractal):
        
            new_subbase = []
            for k in range(self.genus):
        
                index = int(self.base_maps[r][k][0])
                bms = [base[index][0],base[index][-1]]
        
                for m in range(2):
                    if self.base_maps[r][k][-1] == '1':
                        if bms[m][-1] == '1':
                            bms[m] = bms[m][:-1]
                        else:
                            bms[m] = bms[m] + '1'
        
                if self.base_maps[r][k][-1] == '1':
                    new_subbase = new_subbase + list(reversed(bms))
                else:
                    new_subbase = new_subbase + bms
                
            new_base = new_base + [new_subbase]
        
        return new_base

    def get_jun(self,id_bm,n):
        '''Функция находит все стыки на n-ом подразделении'''
        base = self.base_maps

        first = [[[base[r][k]] for k in range(self.genus-1)] for r in range(self.fractal)]
        second = [[[base[r][k+1]] for k in range(self.genus-1)] for r in range(self.fractal)]

        for l in range(n-1):
            base = self.get_bms(base)
            first = [[first[r][k] + [base[r][2*k+1]] 
                      for k in range(self.genus-1)] for r in range(self.fractal)]
            second = [[second[r][k] + [base[r][2*k+2]] 
                       for k in range(self.genus-1)] for r in range(self.fractal)]
        
        for r in range(self.fractal):
            for l in range(self.genus-1):
                for k in range(n-2,-1,-1):
                    first[r][l][-1] = self.get_bm_bm(id_bm,first[r][l][-1],first[r][l][k])
                    second[r][l][-1] = self.get_bm_bm(id_bm,second[r][l][-1],second[r][l][k])

        jun = []
        for r in range(self.fractal):    
            for k in range(self.genus-1):
                jun = jun + [[first[r][k][-1],second[r][k][-1]]]

        jun = set([self.get_con_junction(id_bm,k) for k in jun])

        return jun

    def get_junctions(self):
        '''Функция находит все стыки кривой'''
        id_bm = 'n' + self.alph
        
        # Находим все стыки на первом подразделении
        all_junctions_sub_1 = [[self.base_maps[r][k],self.base_maps[r][k+1]] 
                                for k in range(self.genus-1) for r in range(self.fractal)]
        
        # Приводим стыки на первом подразделении к каноническому виду
        junctions = set([self.get_con_junction(id_bm,k) for k in all_junctions_sub_1])
        
        for k in range(2,20):
            jun = self.get_jun(id_bm,k)
            
            N = len(junctions)
            junctions = junctions.union(jun)
            
            if N == len(junctions):
                break
            
        print('глубина кривой -',k-1)
            
        return sorted(junctions)
    
    def get_all_vertices(self):
        
        proto = self.get_proto()
        
        all_vertices = []
        max_coord = self.div-1
        for l in range(self.fractal):
            vertices = []
            for k in range(self.genus):
                new_vertice = tuple([k for k in proto[l][k] if k in [0,max_coord]])
                if len(new_vertice)==self.dim: vertices.append(new_vertice)
            all_vertices.append(vertices)
        
        return all_vertices
    
    def get_vertex_moments(self):
        
        proto = self.get_proto()            
        all_vertices = self.get_all_vertices()
        
        # Находим номера фракций, которые содержат вершины прототипа
        fractions_numb = []
        for l in range(self.fractal):
            fraction_numb = []
            for k in range(2**self.dim):
                fraction_numb.append(proto[l].index(all_vertices[l][k]))
            fractions_numb.append(fraction_numb)
        
        # Находим номера вершин в подфракциях, попадающих в вершины прототипа
        sub_1 = self.get_subdiv(1,plot=False)
        sub_1_coord = [self.get_curve_coord(k) for k in sub_1]
        subfractions_numb = []
        for l in range(self.fractal):
            subfraction_numb = []
            shift_coord = [min([row[k] for row in sub_1_coord[l]]) for k in range(self.dim)]
            for k in range(2**self.dim):            
                norm_coord = tuple([all_vertices[l][k][m]*self.div/(self.div-1) + shift_coord[m] for m in range(self.dim)])
                D = sub_1_coord[l].index(norm_coord)%(2**self.dim-1)
                subfraction_numb.append(D)
                
            if subfraction_numb[-1]==0: subfraction_numb[-1] = 2**self.dim-1

            subfractions_numb.append(subfraction_numb)
        
        # Составляем систему линейных уравнений и решаем ее
        all_vertex_moments = []
        for l in range(self.fractal):
            A = sp.eye(2**self.dim)
            B = sp.zeros(1,2**self.dim)
            for k in range(2**self.dim):
                A[k,subfractions_numb[l][k]] -= sp.Rational(1,self.genus)
                B[k] = sp.Rational(fractions_numb[l][k],self.genus)
            vertex_moments = A.inv().dot(B)
            all_vertex_moments.append(vertex_moments)

        def lcm(a, b):
            '''Функция находит общий знаменатель для пары чисел'''
            return int(a * b / gcd(a, b))

        def lcms(*numbers):
            '''Функция находит общий знаменатель для списка чисел'''
            return reduce(lcm, numbers)
        
        # Находим знаменатели всех моментов
        all_denominators = []
        for l in range(self.fractal):    
            denominators = [all_vertex_moments[l][k].q for k in range(2**self.dim)]
            all_denominators.append(denominators)

        # Находим общий знаменатель для всех моментов
        all_denominators = sum(all_denominators,[])
        common_denominators = lcms(*all_denominators)
                
        # Приводим моменты к целым числам  
        all_vertex_moments = [[all_vertex_moments[l][k]*common_denominators 
                               for k in range(2**self.dim)] for l in range(self.fractal)]
        
        # Добавляем общий знаменатель к списку моментов
        all_vertex_moments.append(common_denominators)
        
        return all_vertex_moments
