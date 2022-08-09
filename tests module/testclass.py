import itertools as it
import sympy as sp
from functools import reduce
from math import gcd
from numba import njit, float64, int32, prange, int8, boolean, typeof
import numpy as np
import time


class TestClass:
    
    '''
    Основные функции класса:
    is_junction_bad(bm1,bm2) - проверяет является ли стык плохим (для частично заданной кривой), bm1 = 0jIk~
    is_junctions_good(base_maps) - проверяет отсутствие плохих произодных стыков (для полностью заданной кривой), base_maps = [[bm1,bm2,...]]
    
    get_good_curve() - получает список базовых преобразований для кривых без плохих стыков
    
    base_maps() - находит по прототипу базовые преобразования для однозначно определенной кривой
    
    chain_proto = ['jiJijjIIjjiJijiiJIJiJIJi']
    
    В векторной системе
    TestClass(chain_proto,'->ijk').get_base_maps(['0jI~','0ij','0Ji~','0IJ'])
    
    В базисной форме записи
    TestClass(chain_proto,'ijk->').get_base_maps(['0Ji~','0ij','0jI~','0IJ'])
    
    '''
    
    def __init__(self,chain_proto,coding_system,base_maps=None,
                 cut_chain_proto=None,div=None,dim=None,alph=None,
                 fractal=None,genus=None,vect_dict=None):
        
        '''
        chain_proto: указать прототипы в виде списка ['jiJ','jiJ']
        coding_system: указать систему кодирования базовых преобразований ->ijk или ijk->        
        '''
        
        self.chain_proto = chain_proto
        self.coding_system = coding_system
        self.base_maps = None 
        self.cut_chain_proto = self.cut_chain_proto() 
        self.div = div if div is not None else self.get_div()
        self.dim = dim if dim is not None else self.get_dim()
        self.alph = alph if alph is not None else self.get_alph()
        self.fractal = fractal if fractal is not None else self.get_fractal()
        self.genus = genus if genus is not None else self.get_genus()
        self.vect_dict = vect_dict if vect_dict is not None else self.get_vect_dict()
        
    def cut_chain_proto(self):
        '''get prototypes from first and last vectors'''
        return [l[0] + l[-1] for l in self.chain_proto] 
        
    def get_dim(self):
        '''get the curve dimension'''
        return len(set(''.join(self.chain_proto[0]).lower()))
        
    def get_alph(self):
        '''get alphabet for curve prototypes and base maps'''
        return 'ijklmnop'[:self.dim]
    
    def get_fractal(self):
        '''get the curve fractality'''
        return len(self.chain_proto)
    
    def get_div(self):
        '''get the curve div'''
        chain = chain_proto[0]
        if chain != 1:
            chain = ''.join(chain)

        div = 1 
        for k in range(len(chain)):
            if chain[k] == 'i':
                div+=1
            elif chain[k] == 'I':
                div-=1
        return div
    
    def get_genus(self):
        '''get the curve genus'''
        return len(self.chain_proto[0])+1 
    
    def get_vector_code(self,bm):
        '''get vector kode from basis kode of base map'''
        cut_bm = bm[1:self.dim+1]
        
        #Создаем словарь для перекодировки
        dict_bm={}
        for m in range(self.dim):
            dict_bm[self.alph[m]] = cut_bm[m]
            dict_bm[self.alph[m].upper()] = cut_bm[m].swapcase()
            
        inv_dict_bm = {str(m1): k1 for k1, m1 in dict_bm.items()}
        
        #Меняем кодировку (применяем словарь)
        new_bm = ''
        for m in range(self.dim):
            new_bm = new_bm + inv_dict_bm[self.alph[m]]
        
        #Добавляем номер прототипа
        new_bm = bm[0] + new_bm
        
        #Учитываем обращение по времени
        if bm[-1] == '~':
            new_bm = new_bm + '~'
            
        return new_bm
    
    def get_vector_codes(self,base_maps):
        '''get vector kodes from basis kodes of base maps'''
        for l in range(self.fraction):
            for k in range(self.genus):
                base_maps[l][k] = self.get_vector_code(base_maps[l][k])
        
        return base_maps
    
    def get_fraction(self,sub,bm):
        '''apply base map and reverse in vector form to some curve fraction'''
        # Проверяем наличие обращения по времени в базовом преобразовании
        if bm[-1] == '~':
            # Меняем напраления векторов
            bm = bm[:-1].swapcase()
            # Проходим вектора в обратном порядке
            sub = reversed(sub)
            
        # Создаем словарь базового преобразования, например kIJ = {k->i,I->j,J->k} 
        # и его инверсию {K->I,i->j,j->K}
        dict_bm={}
        for k in range(self.dim):
            dict_bm[bm[k]] = self.alph[k]
            dict_bm[bm[k].swapcase()] = self.alph[k].upper()
            
        # Поворачиваем фракцию (применяем словарь)
        fraction = [''.join(list(map(dict_bm.get, k))) for k in sub]
        
        return fraction
    
# Находим все стыки в векторной форме

    def get_con_junction(self,junction):
        '''Функция приводит стык к каноническому виду'''
        bm2 = ''
        id_bm = self.alph
        for l in range(self.dim):
            m = junction[0].lower().index(id_bm[l])
            if junction[0][m] == junction[0][m].upper():
                bm2 = bm2 + junction[1][m].swapcase()
            else:
                bm2 = bm2 + junction[1][m]
                
        if junction[0][-1] == '~': id_bm = id_bm + '~'
        if junction[1][-1] == '~': bm2 = bm2 + '~'
            
        bm1 = junction[0][0] + id_bm
        bm2 = junction[1][0] + bm2
            
        con_junction = (bm1,bm2)
        
        return con_junction

    def get_bm_bm(self,bm2,bm):
        '''Возвращает базовое преобразование базового преобразования. Например bm2(bm) = KiJ(jkI)'''
        new_bm = ''
        id_bm = 'n' + self.alph
        for k in range(1,self.dim+1):
            m = id_bm.index(bm2[k].lower())   
            if bm[m]  == bm[m].upper():
                new_bm = new_bm + bm[m].upper() if bm2[k] == bm2[k].lower() else new_bm + bm[m].lower()
            else:
                new_bm = new_bm + bm[m].upper() if bm2[k] == bm2[k].upper() else new_bm + bm[m].lower() 
        
        if bm[-1] == '~': new_bm = new_bm + '~'
        
        #if bm[-1] == '~' or bm2[-1] == '~':
        #    new_bm = new_bm + '~'
            
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
                    if self.base_maps[r][k][-1] == '~':
                        if bms[m][-1] == '~':
                            bms[m] = bms[m][:-1]
                        else:
                            bms[m] = bms[m] + '~'
        
                if self.base_maps[r][k][-1] == '~':
                    new_subbase = new_subbase + list(reversed(bms))
                else:
                    new_subbase = new_subbase + bms
                
            new_base = new_base + [new_subbase]
        
        return new_base

    def get_time_norm(self,jun):
        '''Функция выполняет нормировку стыков с обращением по времени'''
        if jun[0][-1]=='~' and jun[1][-1]=='~':
            
            jun = [jun[1][:-1],jun[0][:-1]]
        
        elif int(jun[0][0]) > int(jun[1][0]):
            
            if jun[0][-1]=='~' and jun[1][-1]!='~':
                
                jun = [jun[1]+'~',jun[0][:-1]]
                
            elif jun[0][-1]!='~' and jun[1][-1]=='~':
                
                jun = [jun[1][:-1],jun[0]+'~']
            
        return self.get_con_junction(jun)

    def get_jun(self,n):
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
                    first[r][l][-1] = self.get_bm_bm(first[r][l][k],first[r][l][-1])
                    second[r][l][-1] = self.get_bm_bm(second[r][l][k],second[r][l][-1])

        jun = []
        for r in range(self.fractal):    
            for k in range(self.genus-1):
                jun = jun + [[first[r][k][-1],second[r][k][-1]]]

        jun = set([self.get_time_norm(k) for k in jun])

        return jun

    def get_junctions(self):
        '''Функция находит все стыки кривой'''

        # Находим все стыки на первом подразделении
        all_junctions_sub_1 = [[self.base_maps[r][k],self.base_maps[r][k+1]] 
                                for k in range(self.genus-1) for r in range(self.fractal)]
        
        #Нормализуем и приводим стыки к каноническому виду
        junctions = set([self.get_time_norm(k) for k in all_junctions_sub_1])
        
        for k in range(2,20):
            jun = self.get_jun(k)
            
            N = len(junctions)
            junctions = junctions.union(jun)
            
            if N == len(junctions):
                break
            
        #print('глубина кривой -',k-1)
            
        return sorted(junctions)
    
    # Провека на плохие и совершенные стыки
    
    def is_junction_bad(self,bm1,bm2):
        '''check good or bad junction'''
        if self.coding_system == 'ijk->':
            bm1 = self.get_vector_code(bm1)
            bm2 = self.get_vector_code(bm2)
        
        # Применяем бозовое преобразование к cut_chain_proto, где каждый из прототипов
        # состоит из первого и последненго векторов (этого достаточно)
        first_fraction  = self.get_fraction(self.cut_chain_proto[int(bm1[0])],bm1[1:])
        second_fraction = self.get_fraction(self.cut_chain_proto[int(bm2[0])],bm2[1:])
        
        if first_fraction[-1] == second_fraction[0]:
            return True
        else:
            return False
    
    def is_junctions_good(self,base_maps):
        '''Проверяет все стыки кривой'''
        if self.coding_system == 'ijk->':
            base_maps = self.get_vector_codes(base_maps)
        
        self.base_maps = base_maps

        junctions = self.get_junctions()
    
        for m in junctions:    
            if self.is_junction_bad(m[0],m[1]) == True:
                return False
            
        return True        
    
    def is_junction_perfect(self,bm1,bm2):
        '''check perfect junction'''
        if self.coding_system == 'ijk->':
            bm1 = self.get_vector_code(bm1)
            bm2 = self.get_vector_code(bm2)
        
        first_fraction  = self.get_fraction(self.cut_chain_proto[int(bm1[0])],bm1[1:])
        second_fraction = self.get_fraction(self.cut_chain_proto[int(bm2[0])],bm2[1:])
        
        if first_fraction[-1] == second_fraction[0].swapcase():
            return True
        else:
            return False

    def is_junctions_perfect(self,base_maps):
        '''Проверяет все стыки кривой'''
        if self.coding_system == 'ijk->':
            base_maps = self.get_vector_codes(base_maps)
        
        self.base_maps = base_maps

        junctions = self.get_junctions()
    
        for m in junctions:
            if self.is_junction_perfect(m[0],m[1]) == False:
                return False,None,None
            
        return True,junctions,self.get_vertex_moments()

# Находим переходные ломанные, группы вращеиний и хорошие кривые

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

    def get_subdiv(self,sub_numb,plot=True):
        '''get n-th curve subdivision'''
        # Определяем нулевое подразделение кривой
        sub_k = self.chain_proto if plot==True else self.get_vertex_chain_proto()
        
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

    def get_vertex_moments(self):
        '''Функция находит вершинные моменты'''
        proto = self.get_proto()
        all_vertices = self.get_all_vertices()
        
        # Находим номера фракций, которые содержат вершины прототипа
        fractions_numb = []
        for l in range(self.fractal):
            for k in range(2**self.dim):
                fractions_numb.append(proto[l].index(all_vertices[l][k]))
                
        # Находим номера вершин в подфракциях, попадающих в вершины прототипа
        sub_1 = self.get_subdiv(1,plot=False)
        sub_1_coord = [self.get_curve_coord(k) for k in sub_1]
        subfractions_numb = []
        reverse_time = sp.zeros(1,2**self.dim*self.fractal)
        for l in range(self.fractal):
            #Находим сдвиг по всем координатам (для реберных и граневых кривых)
            shift_coord = [min([row[k] for row in sub_1_coord[l]]) for k in range(self.dim)]
            for k in range(2**self.dim):
                #Последовательно сдвигаем и масштабируем каждую вершину прототипа          
                norm_coord = tuple([all_vertices[l][k][m]*self.div/(self.div-1) + 
                                    shift_coord[m] for m in range(self.dim)])
                #Находим номер вершины во фракции, которая сототвествует вершине прототипа
                i = sub_1_coord[l].index(norm_coord)-fractions_numb[k+l*2**self.dim]*(2**self.dim-1)
                #Проверяем есть ли обращение по времени
                if self.base_maps[l][k][-1]=='~':
                    i = 2**self.dim-1 - i
                    reverse_time[k+l*2**self.dim] = 1
                #Определяем к какой кривой принадлежит этот момент
                m = int(self.base_maps[l][k][0])*(2**self.dim)    
                subfractions_numb.append(i+m)

        all_vertex_moments = []
        A = sp.eye(2**self.dim*self.fractal)
        B = sp.zeros(1,2**self.dim*self.fractal)
        for k in range(2**self.dim*self.fractal):
            if reverse_time[k] == 1:
                A[k,subfractions_numb[k]] += sp.Rational(1,self.genus)
                B[k] = sp.Rational(1+fractions_numb[k],self.genus)
            else:    
                A[k,subfractions_numb[k]] -= sp.Rational(1,self.genus)
                B[k] = sp.Rational(fractions_numb[k],self.genus)
        all_vertex_moments = A.inv()*B.T #.dot(B)
        
        def lcm(a, b):
            return int(a * b / gcd(a, b))

        def lcms(*numbers):
            return reduce(lcm, numbers)

        # Находим знаменатели всех моментов
        all_den = [all_vertex_moments[k].q for k in range(2**self.dim*self.fractal)]

        # Находим общий знаменатель для всех моментов
        common_den = lcms(*all_den)

        # Приводим моменты к целым числам  
        all_vertex_moments = [int(all_vertex_moments[k]*common_den) for k in range(2**self.dim*self.fractal)]

        vertex_moments = []
        for l in range(self.fractal):
            vertex_moments.append(all_vertex_moments[l*2**self.dim:(l+1)*2**self.dim])

        vertex_moments.append(common_den)
        
        return vertex_moments
    
    def get_all_vertices(self):
        '''Функция находит координаты вершин прототипов, в порядке их прохождения'''
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
    
    def get_vertex_chain_proto(self):
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
    
    def sum_coord(self,coord1,coord2):
            return tuple(a_n + b_n for a_n, b_n in zip(coord1, coord2))
        
    def get_grid(self,proto):
        '''get grids div**dim from protypies vertices'''
        #Определяем сдвиги координат для разных размерностей
        if self.dim == 2:
            shift = ['i','j','ij']
        elif self.dim == 3:
            shift = ['i','j','k','ij','jk','ik','ijk']
        
        #Строим сетку
        grid = []
        for k in range(self.genus):
    
            coord_frac = []
            coord_frac.append(proto[k])
            for m in range(len(shift)):
                coord_frac.append(self.sum_coord(proto[k],self.vect_dict[shift[m]]))

            grid.append(set(coord_frac))
        
        #Определяем тип сетки
        one_side_coord = tuple([self.div-1]+[0]*(self.dim-1))
        if proto[-1] == one_side_coord:
            type_grid = 'one_side'
            
        diag_coord = tuple([self.div-1]*self.dim)
        if proto[-1] == diag_coord:
            type_grid = 'diag'

        return grid,type_grid

    def get_shift_one_side(self):
        '''Находим сдвиги координат для односторонних кривых'''
        one_side = self.alph + self.alph.upper()
        one_side = [k for k in one_side]
        
        return one_side
    
    def get_shift_diag(self):
        '''Находим сдвиги координат для диагональных кривых'''
        all_letters = [[self.alph[k],self.alph[k].upper()] for k in range(self.dim)]
        all_comb_let = list(it.combinations(all_letters,self.dim))
        diag = list(map(''.join, it.chain(*[it.product(*k) for k in all_comb_let])))
        
        return diag

    def get_tran_broken(self,grid,type_grid,start=None,end=None):
        '''get transitional broken line on the grid''' #Пока только одна ломанная
        #Находим начальную координату
        if start == None:
            start_coord = tuple([0]*self.dim)
        
        #Находим конечную координату
        if end == None:
            if type_grid == 'one_side':
                end_coord = tuple([self.div]+[0]*(self.dim-1))
            elif type_grid == 'diag':
                end_coord = tuple([self.div]*self.dim)        
        
        #Находим сдвиги координат
        if self.fractal==2:
            step = self.get_shift_one_side() + self.get_shift_diag() 
                
        elif type_grid == 'one_side':
            step = self.get_shift_one_side()
            
        elif type_grid == 'diag':
            step = self.get_shift_diag() 
        
        current_branchs = [[start_coord]]
        vect_current_branchs = [[]]
        #Проходим все координаты
        for l in range(self.genus-1):

            new_branchs = []
            vect_new_branchs = []
            #Проходим все ветви
            for k in range(len(current_branchs)):
                
                Start = current_branchs[k][-1]
                for m in range(len(step)):
                    
                    End = self.vect_dict[step[m]]
                    current_coord = self.sum_coord(Start,End)
                    
                    if current_coord in grid[l]:
                        if current_coord in grid[l+1]:
        
                            new_branchs.append(current_branchs[k]+[current_coord])
                            vect_new_branchs.append(vect_current_branchs[k]+[step[m]])
                        
            current_branchs = new_branchs  
            vect_current_branchs = vect_new_branchs
            
        new_branchs = []
        vect_new_branchs = []
        for k in range(len(current_branchs)):
            
            Start = current_branchs[k][-1]
            for m in range(len(step)):
                
                End = self.vect_dict[step[m]]
                current_coord = self.sum_coord(Start,End)
                
                if current_coord == end_coord:
                    
                    new_branchs.append(current_branchs[k]+[end_coord])
                    vect_new_branchs.append(vect_current_branchs[k]+[step[m]])
        
        tran_broken = new_branchs
        vect_tran_broken = vect_new_branchs
        
        return tran_broken,vect_tran_broken

    def get_bms_fraction(self,rev=True):
        '''get list base maps with reversed'''
        all_letters = [[self.alph[k],self.alph[k].upper()] for k in range(self.dim)]
        all_comb_let = list(it.permutations(all_letters,self.dim))
        bms = list(map(''.join, it.chain(*[it.product(*k) for k in all_comb_let])))
        
        if rev == True:    
            bms = bms + [k+'~' for k in bms]
        
        return bms

    def get_rotation_groups(self,proto,rev=False,bms=None):
        
        grid,type_grid = self.get_grid(proto)    
   
        tran_broken,vect_tran_broken = self.get_tran_broken(grid,type_grid)
        
        vertex_chain_proto = self.get_vertex_chain_proto()
        
        if bms==None:
            bms = self.get_bms_fraction(rev)
        
        all_groups = []
        for l in range(len(tran_broken)):
        
            groups = []
            for m in range(self.genus):
                
                #Проверяем является кривая односторонней или диагональной
                n = '0' if len(vect_tran_broken[l][m])==1 else '1'
                
                group = []
                for k in range(len(bms)):
                    
                    bm_fraction = self.get_fraction(vertex_chain_proto[int(n)],bms[k])
                    coord_frac = self.get_curve_coord(bm_fraction,tran_broken[l][m])
                    if coord_frac[-1] == tran_broken[l][m+1]:
                        if set(coord_frac) == grid[m]:
                            
                            group.append(n + bms[k])
        
                groups.append(tuple(group))
            
            all_groups.append(groups)
            
        return all_groups
    
    def get_base_maps(self,bms):
        
        if self.coding_system == 'ijk->':
            bms = [self.get_vector_code(k) for k in bms]
        
        proto = self.get_proto()[0]

        groups = self.get_rotation_groups(proto,bms)
        
        base_maps = [j[0] for j in groups]
        
        if self.coding_system == 'ijk->':
            base_maps = [self.get_vector_code(k) for k in base_maps]
        
        return base_maps
        
    def gef_perfect_curve(self,groups):
                
        current_branchs = [[j] for j in groups[0]]

        for l in range(1,self.genus):
    
            new_branchs = []
            for k in range(len(current_branchs)):
    
                for m in range(len(groups[l])):
        
                    if self.is_junction_perfect(current_branchs[k][-1],groups[l][m]) == True:
                        
                        new_branchs.append(current_branchs[k]+[groups[l][m]])
        
            current_branchs = new_branchs

        return current_branchs      
        
    def get_good_curve(self):
        ''''get base maps for good curve'''
        proto = self.get_proto()[0]
        groups = self.get_rotation_groups(proto)
        
        t = time.time()
        
        current_branchs = [[j] for j in groups[0]]

        for l in range(1,self.genus):
    
            new_branchs = []
            for k in range(len(current_branchs)):
    
                for m in range(len(groups[l])):
        
                    if self.is_junction_bad(current_branchs[k][-1],groups[l][m]) == False:
                        
                        new_branchs.append(current_branchs[k]+[groups[l][m]])
        
            current_branchs = new_branchs
        
        elapsed = time.time() - t   
        
        print(len(current_branchs),'- кривые с хорошими стыками на первом подразделении')
        print('время перебора -', elapsed)
        
        t = time.time()
        
        for k in range(len(current_branchs)-1,-1,-1):
    
            if self.is_junctions_good([current_branchs[k]]) == False:
                del current_branchs[k]
                
        elapsed = time.time() - t   
        
        print(len(current_branchs),'- кривые с хорошими производными стыками') 
        print('время перебора -', elapsed)
        
        return current_branchs
    

def get_moment_sub_k(moments,n):
    
    den = moments[-1]
    moment_sub_k = [np.array(moments[k]) for k in range(self.fractal)]

    for n in range(n):

        last_moments = []

        for r in range(self.fractal):
            current_moment = []
    
            for k in range(self.genus):
        
                new_moment = moment_sub_k[int(self.base_maps[r][k][0])]
        
                if self.base_maps[r][k][-1]=='~':
                    new_moment = np.flip(den*self.genus**n - new_moment)
        
                norm_moment = new_moment + k*den*self.genus**n
                current_moment = np.concatenate((current_moment, norm_moment))
            
            last_moments.append(current_moment)
    
        moment_sub_k = last_moments
        
    moment_sub_k.append(den)
        
    return [int32(k) for k in moment_sub_k]



def get_junction_ratio(junction,moments,n):
       
    bm1,bm2 = junction
    
    # Генерируем моменты n-ого подразделения 
    moment_sub_k = get_moment_sub_k(moments,n)
    
    first_moments = moment_sub_k[int(bm1[0])]
    second_moments = moment_sub_k[int(bm2[0])]
    
    if bm1[-1] == '~':
        first_moments = np.flip(moment_sub_k[-1]*self.genus**n - first_moments)
        
    if bm2[-1] == '~':
        second_moments = np.flip(moment_sub_k[-1]*self.genus**n - second_moments)
    
    d = np.arange(0, 2**self.dim*self.genus**(n), 2**self.dim)
    d = np.flip(d)
        
    first_moments = np.delete(first_moments, d, 0)
    second_moments = np.delete(second_moments, d, 0)
    
    if n != 0:    
        second_moments = moment_sub_k[-1]*self.genus**n + second_moments
    else:
        second_moments = moment_sub_k[-1] + second_moments

    #print(first_moments)    

    # Генерируем n-ое подразделение кривой
    sub = self.get_subdiv(n,plot=False)

    frac_1 = self.get_fraction(sub[int(bm1[0])],bm1[1:])
    coord_1 = self.get_curve_coord(frac_1)
    
    frac_2 = self.get_fraction(sub[int(bm2[0])],bm2[1:])
    coord_2 = self.get_curve_coord(frac_2,start=coord_1[-1])
    
    coord_1.pop(0)
    coord_2.pop(0)
    
    coord_1 = np.array(coord_1)
    coord_2 = np.array(coord_2)
    
    x1 = coord_1
    x2 = coord_2

    t1 = first_moments
    t2 = second_moments
    
    r = get_ratio_l1_3d(x1,x2,t1,t2)
        
    return moment_sub_k[-1]*r



def get_matrix_base_maps():
    
    d = { 'i':1, 'I':-1, 'j':2, 'J':-2, 'k':3, 'K':-3, '0':0, '1':1 }

    matrix_base_maps = np.zeros((self.fractal,self.genus,self.dim+2),dtype=np.int8)
    for l in range(self.fractal):
        for k in range(self.genus):
            for m in range(self.dim+1):
                matrix_base_maps[l][k][m] = d[self.base_maps[l][k][m]]

    for l in range(self.fractal):
        for k in range(self.genus):
            if self.base_maps[l][k] == '~':
                matrix_base_maps[l][k][-1] = 1
            else:
                matrix_base_maps[l][k][-1] = 0
    
    return matrix_base_maps

def get_all_matrix_base_maps(base_maps):
    
    d = { 'i':1, 'I':-1, 'j':2, 'J':-2, 'k':3, 'K':-3, '0':0, '1':1 }

    all_matrix_base_maps = np.zeros((len(base_maps),self.genus,self.dim+2),dtype=np.int8)
    for l in range(len(base_maps)):
        for k in range(self.genus):
            for m in range(self.dim+1):
                all_matrix_base_maps[l][k][m] = d[base_maps[l][k][m]]

    for l in range(len(base_maps)):
        for k in range(self.genus):
            if base_maps[l][k] == '~':
                all_matrix_base_maps[l][k][-1] = 1
            else:
                all_matrix_base_maps[l][k][-1] = 0
    
    return all_matrix_base_maps




dim = 3
fractal = 2
genus = 8


@njit(int8[:,:](int8[:,:]),cache=True)
def get_con_junction_new(jun):
    
    bm2 = np.zeros(dim+2,dtype=np.int8)
    id_bm = np.array([9,1,2,3,0],dtype=np.int8)
    
    for l in range(1,dim+1):
        m = np.where(np.abs(jun[0][1:dim+1]) == id_bm[l])[0][0]+1

        if jun[0][m] == -abs(jun[0][m]):
            bm2[l] = -1*jun[1][m]
        else:
            bm2[l] = jun[1][m]
    
    if jun[0][-1] == 1: id_bm[-1] = 1
    if jun[1][-1] == 1: bm2[-1] = 1
            
    id_bm[0] = jun[0][0]
    bm2[0] = jun[1][0]
    
    jun[0] = id_bm
    jun[1] = bm2
    
    return jun

@njit(int8[:,:](int8[:,:]),cache=True)
def get_time_norm_new(jun):
    #Функция выполняет нормировку стыков с обращением по времени
    if jun[0][-1]==1 and jun[1][-1]==1:
        
        new_jun = np.zeros((2,dim+2),dtype=np.int8)
        
        new_jun[0] = jun[1]
        new_jun[1] = jun[0]
        
        new_jun[0][-1] = 0
        new_jun[1][-1] = 0
        
        jun = new_jun
        
    elif jun[0][0] > jun[1][0]:
            
        if jun[0][-1]==1 and jun[1][-1]==0:
            
            new_jun = np.zeros((2,dim+2),dtype=np.int8)
            
            new_jun[0] = jun[1]
            new_jun[1] = jun[0]
            
            new_jun[0][-1] = 1
            new_jun[1][-1] = 0
        
            jun = new_jun
                
        elif jun[0][-1]==0 and jun[1][-1]==1:
            
            new_jun = np.zeros((2,dim+2),dtype=np.int8)
                
            new_jun[0] = jun[1]
            new_jun[1] = jun[0]
            
            new_jun[0][-1] = 0
            new_jun[1][-1] = 1
        
            jun = new_jun
            
    return get_con_junction_new(jun)

@njit(int8[:](int8[:],int8[:]),cache=True)
def get_bm_bm_new(bm2,bm):

    new_bm = np.zeros(dim+2,dtype=np.int8)
    id_bm = np.array([9,1,2,3],dtype=np.int8)
    for k in range(1,dim+1):
        m = np.where(id_bm == abs(bm2[k]))[0][0]
        if bm[m]  == -abs(bm[m]):
            new_bm[k] = -abs(bm[m]) if bm2[k] == abs(bm2[k]) else abs(bm[m])
        else:
            new_bm[k] = -abs(bm[m]) if bm2[k] == -abs(bm2[k]) else (bm[m])

        if bm[-1] == 1:
            new_bm[-1] = 1

        new_bm[0] = bm[0]
        
    return new_bm

@njit(int8[:,:,:](int8[:,:,:],int8[:,:,:]),cache=True)
def get_bms_new(base,matrix_base_maps):

    new_base = np.zeros((fractal,2*genus,dim+2),dtype=np.int8)
    for r in range(fractal):
        for k in range(genus):
        
            index = matrix_base_maps[r][k][0]
            bms = np.zeros((2,dim+2),dtype=np.int8)
            bms[0] = base[index][0]
            bms[1] = base[index][-1]

            for m in range(2):
                if matrix_base_maps[r][k][-1] == 1:
                    if bms[m][-1] == 1:
                        bms[m][-1] = 0
                    else:
                        bms[m][-1] = 1

            if matrix_base_maps[r][k][-1] == 1:
                new_base[r][2*k] = bms[1]
                new_base[r][2*k+1] = bms[0]
            else:
                new_base[r][2*k] = bms[0]
                new_base[r][2*k+1] = bms[1]
                
    return new_base

@njit(int8[:,:,:](int8[:,:,:],int8),cache=True)
def get_jun_new(matrix_base_maps,n):
    
    base = matrix_base_maps
    
    first = np.zeros((fractal,genus-1,n,dim+2),dtype=np.int8)
    second = np.zeros((fractal,genus-1,n,dim+2),dtype=np.int8)
    for r in range(fractal):
        for k in range(genus-1):
            first[r][k][0] = base[r][k]
            second[r][k][0] = base[r][k+1]
    
    for l in range(n-1):
        base = get_bms_new(base,matrix_base_maps)
        for r in range(fractal):
            for k in range(genus-1):
                first[r][k][l+1] = base[r][2*k+1]
                second[r][k][l+1] = base[r][2*k+2]

    for r in range(fractal):
        for l in range(genus-1):
            for k in range(n-2,-1,-1):
                first[r][l][-1] = get_bm_bm_new(first[r][l][k],first[r][l][-1])
                second[r][l][-1] = get_bm_bm_new(second[r][l][k],second[r][l][-1])
    
    jun = np.zeros((fractal*(genus-1),2,dim+2),dtype=np.int8)
    for r in range(fractal):    
        for k in range(genus-1):
            jun[(genus-1)*r+k][0] = first[r][k][-1]
            jun[(genus-1)*r+k][1] = second[r][k][-1]
            jun[(genus-1)*r+k] = get_time_norm_new(jun[(genus-1)*r+k])
    
    return jun

@njit(int8[:,:,:](int8[:,:,:]),cache=True)
def get_set(jun):
    
    for m in range(len(jun)-1,-1,-1):
        for k in range(m-1,-1,-1):
            if np.array_equal(jun[m],jun[k]):
                
                #Не работает третий аргумент delete
                #jun = np.delete(jun, m, int8(0))
                #Удаление стыка с помощью логической маски
                mask = np.zeros(jun.shape[0], dtype=np.int8) == 0
                mask[m] = False
                jun = jun[mask]
       
                break
    
    return jun

@njit(int8[:,:,:](int8[:,:,:]),cache=True)
def get_junctions_new(matrix_base_maps):
    
    junctions = np.zeros((fractal*(genus-1),2,dim+2),dtype=np.int8)
    for r in range(fractal):
        for k in range(genus-1):
            junctions[(genus-1)*r+k][0] = matrix_base_maps[r][k]
            junctions[(genus-1)*r+k][1] = matrix_base_maps[r][k+1]

    for k in range(fractal*(genus-1)):
        junctions[k] = get_time_norm_new(junctions[k]) 

    junctions = get_set(junctions)

    for k in range(2,20):
        jun = get_jun_new(matrix_base_maps,k)

        N = len(junctions)
        junctions = np.concatenate((junctions, jun), axis=0)
    
        junctions = get_set(junctions)

        if N == len(junctions):
            break
    
    #print('глубина кривой -',k-1)

    return junctions

def get_conv_jun(jun,dim):

    d = { 1:'i', -1:'I', 2:'j', -2:'J', 3:'k', -3:'K'}

    conv_jun = []
    for k in jun:
        bm1 = str(k[0][0]) + 'ijk'
        
        if str(k[0][-1]) == 1:
            bm1 = bm1 + '~'         
        
        bm2 = str(k[1][0]) + ''.join(map(d.get, k[1][1:dim+1]))
    
        if str(k[1][-1]) == 1:
            bm2 = bm2 + '~' 

        conv_jun.append(tuple([bm1,bm2]))
        
    return conv_jun

def is_junctions_perfect_new(base_maps):
    #Проверяет все стыки кривой#
    if self.coding_system == 'ijk->':
        base_maps = self.get_vector_codes(base_maps)
        
    self.base_maps = base_maps

    matrix_base_maps = get_matrix_base_maps()
    junctions = get_junctions_new(matrix_base_maps)
    junctions = get_conv_jun(junctions,dim)
    
    for m in junctions:
        if self.is_junction_perfect(m[0],m[1]) == False:
            return False,None,None
            
    return True,junctions,self.get_vertex_moments()


chain_proto = ['kjKikJK','kiKjIki'] #(Токарев-Хаверкорт)
#chain_proto = ['kijIKiJ','kiKjIki'] #(Шалыга-Хаверкорт)

Test = TestClass(chain_proto,'->ijk')
self = Test



@njit(int8[:,:](int8[:,:],int8[:]),cache=True)
def get_fraction_new(sub,bm):

    if bm[-1] == 1:
        # Меняем напраления векторов
        bm = -bm[:-1]
        # Проходим вектора в обратном порядке
        sub = sub[::-1,:]

    # Поворачиваем фракцию
    sub = np.column_stack((np.sign(bm[0])*sub[:,abs(bm[0])-1],
                           np.sign(bm[1])*sub[:,abs(bm[1])-1],
                           np.sign(bm[2])*sub[:,abs(bm[2])-1])) 
    
    return sub

@njit(boolean(int8[:,:,:]),cache=True)
def is_junction_perfect_new(junctions):
    
    #Кривая Токорева
    cut_chain_proto_np = np.array([[[0,0,1],[0,0,-1]],
                                   [[0,0,1],[1,0,0]]],dtype=np.int8)
    
    #Кривая Шалыги
    #cut_chain_proto_np = np.array([[[0,0,1],[0,-1,0]],
    #                               [[0,0,1],[1,0,0]]],dtype=np.int8)
    
    for k in range(len(junctions)):
        
        bm1 = junctions[k][0]
        bm2 = junctions[k][1]
        
        first_fraction  = get_fraction_new(cut_chain_proto_np[bm1[0]],bm1[1:])
        second_fraction = get_fraction_new(cut_chain_proto_np[bm2[0]],bm2[1:])
        
        if np.array_equal(first_fraction[-1],-second_fraction[0])==False:
            return False
        
    return True


proto = self.get_proto()[0]
grid,type_grid = self.get_grid(proto)
tran_broken,vect_tran_broken = self.get_tran_broken(grid,type_grid)
rotation_groups_1 = self.get_rotation_groups(proto)#,True

curve_1 = []
for k in rotation_groups_1: 
    curve_1 += self.gef_perfect_curve(k)
print(len(curve_1),'- кривые с совершенными стыками на первом подразделении')

proto = self.get_proto()[1]
grid,type_grid = self.get_grid(proto)
tran_broken,vect_tran_broken = self.get_tran_broken(grid,type_grid)
rotation_groups_2 = self.get_rotation_groups(proto)#,True

curve_2 = []
for k in rotation_groups_2: 
    curve_2 += self.gef_perfect_curve(k)
print(len(curve_2),'- кривые с совершенными стыками на первом подразделении')


x = 0
for k in range(len(rotation_groups_1)):
    g = 1
    for m in range(self.genus):
        g *= len(rotation_groups_1[k][m])
    x += g 

y = 0
for k in range(len(rotation_groups_2)):
    g = 1
    for m in range(self.genus):
        g *= len(rotation_groups_2[k][m])
    y += g 

print('общее количество кривых - ',x*y)




def search_curve(curve_1,curve_2):
    
    c1 = get_all_matrix_base_maps(curve_1)
    c2 = get_all_matrix_base_maps(curve_2)

    all_base_maps = []
    all_junctions = []
    for k in range(len(c1)):
        #t = time.time()
        for m in range(len(c2)):

            junctions = get_junctions_new(np.array([c1[k],c2[m]],dtype=np.int8))

            if is_junction_perfect_new(junctions)==True:
                all_junctions.append(get_conv_jun(junctions,dim))
                all_base_maps.append([curve_1[k],curve_2[m]])
                
        #elapsed = time.time() - t
        #print(k,' - ',elapsed)

    return all_base_maps,all_junctions


t = time.time()

all_base_maps,all_junctions = search_curve(curve_1,curve_2)

elapsed = time.time() - t
print('время перебора -', elapsed)


all_vertex_moments = []
for k in range(len(all_base_maps)):
    self.base_maps = all_base_maps[k]
    all_vertex_moments.append(self.get_vertex_moments())
    


sub_numb = 2
n = (2**self.dim-1)*self.genus**sub_numb

@njit(float64(int32[:,:],int32[:,:],int32[:],int32[:]),parallel=True,cache=True)
def get_ratio_l1_3d(x1,x2,t1,t2):
        
    z = np.zeros(n)
    for m in prange(n):
    
        z[m] = np.max((np.absolute(x1[m,0]-x2[:,0])+
                       np.absolute(x1[m,1]-x2[:,1])+
                       np.absolute(x1[m,2]-x2[:,2]))**3/
                       (t2-t1[m]));
            
    return np.max(z)


@njit(float64(int32[:,:],int32[:,:],int32[:],int32[:]),parallel=True,cache=True)
def get_ratio_l2_3d(x1,x2,t1,t2):
        
    z = np.zeros(n)
    for m in prange(n):
    
        z[m] = np.max(((x1[m,0]-x2[:,0])**2+
                       (x1[m,1]-x2[:,1])**2+
                       (x1[m,2]-x2[:,2])**2)**(3/2)/
                       (t2-t1[m]));
            
    return np.max(z)


t = time.time()

all_ratio = []
for m in range(len(all_base_maps)):
    
    self.base_maps = all_base_maps[m]

    ratio = [get_junction_ratio(k,all_vertex_moments[m],sub_numb) for k in all_junctions[m]]
    
    all_ratio.append(max(ratio))


elapsed = time.time() - t
print('время вычисления l1 - ',elapsed)




