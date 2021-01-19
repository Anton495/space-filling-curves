import itertools as it

class TestClass:
    
    '''
    Основные функции класса:
    check_junction(bm1,bm2) - проверяет хороший или плохой стык (для частично заданной кривой)
    check_all_junctions(base_maps) - проверяет все стыки (для полностью заданной кривой)
    '''
    
    def __init__(self,chain_proto,coding_system,base_maps=None,
                 cut_chain_proto=None,div=None,dim=None,alph=None,
                 fractal=None,genus=None):
        
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

    def check_junction(self,bm1,bm2):
        '''check good or bad junction'''
        if self.coding_system == 'ijk->':
            bm1 = self.get_vector_code(bm1)
            bm2 = self.get_vector_code(bm2)
        
        # Применяем бозовое преобразование к cut_chain_proto, где каждый из прототипов
        # состоит из первого и последненго векторов (этого достаточно)
        first_fraction  = self.get_fraction(self.cut_chain_proto[int(bm1[0])],bm1[1:])
        second_fraction = self.get_fraction(self.cut_chain_proto[int(bm2[0])],bm2[1:])
        
        if first_fraction[-1] == second_fraction[0]:
            return False
        else:
            return True
    
#Находим все стыки в векторной форме

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

    def check_all_junctions(self,base_maps):
        '''Проверяет все стыки кривой'''
        if self.coding_system == 'ijk->':
            base_maps = self.get_vector_codes(base_maps)
        
        self.base_maps = base_maps

        junctions = self.get_junctions()
    
        for m in junctions:    
            if self.check_junction(m[0],m[1]) == False:
                return False
            
        return True
