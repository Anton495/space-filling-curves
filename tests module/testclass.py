class TestClass:
    
    '''
    Основные функции класса:
    check_junction(bm1,bm2) - определяет хороший или плохой стык
    '''
    
    def __init__(self,coding_system,chain_proto,
                 cut_chain_proto=None,dim=None,alph=None,fractal=None,genus=None):
        
        '''
        coding_system: указать систему кодирования базовых преобразований ->ijk или ijk->
        chain_proto: указать прототипы в виде списка ['jiJ','jiJ'] (без диагональных шагов)
        div: адичность кривой
        '''
        
        self.coding_system = coding_system
        self.chain_proto = chain_proto
        self.cut_chain_proto = self.cut_chain_proto()        
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
    
    def get_genus(self):
        '''get the curve genus'''
        return len(self.chain_proto[0])+1 
    
    def get_basis_coding(self,bm):
        '''get basis base map from vector base map'''
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
        
        #Добавляем номер кривой
        new_bm = bm[0] + new_bm
        
        #Учитываем обращение по времени
        if bm[-1] == '~':
            new_bm = new_bm + '~'
            
        return new_bm
    
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
            bm1 = self.get_basis_coding(bm1)
            bm2 = self.get_basis_coding(bm2)
            
        first_fraction  = self.get_fraction(self.cut_chain_proto[int(bm1[0])],bm1[1:])
        second_fraction = self.get_fraction(self.cut_chain_proto[int(bm2[0])],bm2[1:])
        
        if first_fraction[-1] == second_fraction[0]:
            return False
        else:
            return True
    
    def check_base_maps(self,base_maps):
        '''check good or bad all junctions of subdivision 1'''
        for l in range(self.fractal):
            for k in range(self.genus-1):
                if self.check_junction(base_maps[l][k],base_maps[l][k+1]) == False:
                    return False
        
        return True
    
    


    


#example 1

#chain_proto = ['jiJ']
#base_maps = [['0ji','0ij','0ij','0JI']]

#example 2

chain_proto = ['jjiJJijj']
base_maps = [['0ij','0Ij','0ij','0iJ','0IJ','0iJ','0ij','0Ij','0ij']]

Test = TestClass('->ijk',chain_proto)
print(Test.check_base_maps(base_maps))
