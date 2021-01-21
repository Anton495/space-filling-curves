import itertools as it
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
        
    def is_junction_perfect(self,bm1,bm2):
        '''check good or bad junction'''
        if self.coding_system == 'ijk->':
            bm1 = self.get_vector_code(bm1)
            bm2 = self.get_vector_code(bm2)
        
        first_fraction  = self.get_fraction(self.cut_chain_proto[int(bm1[0])],bm1[1:])
        second_fraction = self.get_fraction(self.cut_chain_proto[int(bm2[0])],bm2[1:])
        
        if first_fraction[-1] == second_fraction[0].swapcase():
            return True
        else:
            return False
    
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

        jun = set([self.get_con_junction(k) for k in jun])

        return jun

    def get_time_norm(self,jun0):
        '''Функция выполняет нормировку стыков с обращением по времени'''
        if jun0[0][-1]=='~' and jun0[1][-1]=='~':
            
            jun1 = list(reversed(jun0))
            jun2 = [jun1[0][:-1],jun1[1][:-1]]
            
            new_jun = self.get_con_junction(jun2)
            
            return new_jun
        
        elif int(jun0[0][0]) > int(jun0[1][0]):
            
            if jun0[0][-1]=='~' and jun0[1][-1]!='~':
                
                jun1 = list(reversed(jun0))
                jun2 = [jun1[0],jun1[1][:-1]]
                jun3 = [jun2[0]+'~',jun2[1]]
                
                new_jun = self.get_con_junction(jun3)
                
                return new_jun
                
            elif jun0[0][-1]!='~' and jun0[1][-1]=='~':
                
                jun1 = list(reversed(jun0))
                jun2 = [jun1[0][:-1],jun1[1]]
                jun3 = [jun2[0],jun2[1]+'~']
                
                new_jun = self.get_con_junction(jun3)
                
                return new_jun
            
        return jun0

    def get_junctions(self):
        '''Функция находит все стыки кривой'''

        # Находим все стыки на первом подразделении
        all_junctions_sub_1 = [[self.base_maps[r][k],self.base_maps[r][k+1]] 
                                for k in range(self.genus-1) for r in range(self.fractal)]
        
        # Приводим стыки на первом подразделении к каноническому виду
        junctions = set([self.get_con_junction(k) for k in all_junctions_sub_1])
        # Нормализуем стыки
        junctions = set([self.get_time_norm(k) for k in junctions])

        
        for k in range(2,20):
            jun = self.get_jun(k)
            
            N = len(junctions)
            junctions = junctions.union(jun)
            
            # Нормализуем стыки
            junctions = set([self.get_time_norm(k) for k in junctions])
            
            if N == len(junctions):
                break
            
        #print('глубина кривой -',k-1)
            
        return sorted(junctions)

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

    def get_tran_broken(self,grid,type_grid,start=None,end=None):
        '''get transitional broken line on the grid''' #Пока только одна ломанная
        if start == None:
            start_coord = tuple([0]*self.dim) 

        if type_grid == 'one_side':
            step = ['i','j','I','J'] 
            if end == None:
                end_coord = tuple([self.div]+[0]*(self.dim-1))

        elif type_grid == 'diag':
            step = ['ij','iJ','Ij','IJ']
            if end == None:
                end_coord = tuple([self.div]*self.dim)
    
        tran_broken = [start_coord]
        for k in range(self.genus-1):
    
            start = tran_broken[k]
            for m in range(len(step)):
        
                end = self.vect_dict[step[m]]
                current_coord = self.sum_coord(start,end)
        
                if current_coord in grid[k]:
                    if current_coord in grid[k+1]:
                        tran_broken.append(current_coord)
                        break

        tran_broken.append(end_coord)
        
        return tran_broken

    def get_rotation_groups(self,proto,bms=None):
        
        grid,type_grid = self.get_grid(proto)    
   
        tran_broken = self.get_tran_broken(grid,type_grid)
        
        vertex_chain_proto = self.get_vertex_chain_proto()
        
        if bms==None:
            bms = self.get_bms_fraction()

        groups = []
        for m in range(self.genus):
    
            group = []
            for k in range(len(bms)):
            
                bm_fraction = self.get_fraction(vertex_chain_proto[int(bms[k][0])],bms[k][1:])
                coord_frac = self.get_curve_coord(bm_fraction,tran_broken[m])
                if coord_frac[-1] == tran_broken[m+1]:
                    if set(coord_frac) == grid[m]:
                        group.append(bms[k])
        
            groups.append(tuple(group))
            
        return groups
        
    def get_bms_fraction(self,rev=True):
        '''get list base maps with reversed'''
        all_letters = [[self.alph[k],self.alph[k].upper()] for k in range(self.dim)]
        all_comb_let = list(it.permutations(all_letters,self.dim))
        bms = list(map(''.join, it.chain(*[it.product(*k) for k in all_comb_let])))
        
        bms = ['0'+k for k in bms]
        
        if rev == True:    
            bms = bms + [k+'~' for k in bms]
        
        return bms
    
    def get_base_maps(self,bms):
        
        if self.coding_system == 'ijk->':
            bms = [self.get_vector_code(k) for k in bms]
        
        proto = self.get_proto()[0]

        groups = self.get_rotation_groups(proto,bms)
        
        base_maps = [j[0] for j in groups]
        
        if self.coding_system == 'ijk->':
            base_maps = [self.get_vector_code(k) for k in base_maps]
        
        return base_maps
        
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
