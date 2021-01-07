def convert_bm(base_maps):
    '''Эта функция переводит базовые пробразования вида ijk -> Jki~ в Jki~ -> ijk и наоборот'''
    
    alph = 'ijklmnop'

    fractal = len(base_maps)
    
    genus = len(base_maps[0])
    
    dim = len(base_maps[0][0])
    if fractal != 1:
        dim -= 1
    if base_maps[0][0][-1]=='~':
        dim -= 1
    
    if fractal != 1: 
        shift = 1
    else:
        shift = 0
    
    new_base_maps = []
    for l in range(fractal):
        new_bms = []
        for k in range(genus):
            
            cut_bm = base_maps[l][k][shift:dim+shift]

            dict_bm={}
            for m in range(dim):
                dict_bm[alph[m]] = cut_bm[m]
                dict_bm[alph[m].upper()] = cut_bm[m].swapcase()

            inv_dict_bm = {str(m1): k1 for k1, m1 in dict_bm.items()}

            new_bm = ''
            for m in range(dim):
                new_bm = new_bm + inv_dict_bm[alph[m]]

            new_bm = base_maps[l][k][0] + new_bm

            if base_maps[l][k][-1] == '~':
                new_bm = new_bm + '1'
                
            if base_maps[l][k][-1] == '1':
                new_bm = new_bm + '~'
    
            new_bms.append(new_bm)
        new_base_maps.append(new_bms)
    
    return new_base_maps
