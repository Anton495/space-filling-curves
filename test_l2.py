from examples import *
import numpy as np
import time
from numba import njit, float64, int32, prange

'''Тестовый модуль для вычисления стыков в l2'''

#curve = get_hilbert_curve()
#curve = get_beta_Omega_curve()
#curve = get_ARW_curve()

#curve = get_haverkort_curve_1()
#curve = get_haverkort_curve_2()
#curve = get_tokarev_curve()
#curve = get_neptunus_curve()
#curve = get_luna_curve()

curve = get_17_curve()




sub_numb = 2



n = (2**curve.dim-1)*curve.genus**sub_numb

@njit(float64(int32[:,:],int32[:,:],int32[:],int32[:]),parallel=True,cache=True)
def get_ratio_l2_3d(x1,x2,t1,t2):
        
    z = np.zeros(n)
    for m in prange(n):
    
        z[m] = np.max(((x1[m,0]-x2[:,0])**2+
                       (x1[m,1]-x2[:,1])**2+
                       (x1[m,2]-x2[:,2])**2)**(3/2)/
                       (t2-t1[m]));
            
    return np.max(z)


@njit(float64(int32[:,:],int32[:,:],int32[:],int32[:]),parallel=True,cache=True)
def get_ratio_l2_2d(x1,x2,t1,t2):
        
    z = np.zeros(n)
    for m in prange(n):
    
        z[m] = np.max(((x1[m,0]-x2[:,0])**2+
                       (x1[m,1]-x2[:,1])**2)/
                       (t2-t1[m]));

    return np.max(z)


def get_moment_sub_k(moments,n):
    
    den = moments[-1]
    moment_sub_k = [np.array(moments[k]) for k in range (curve.fractal)]

    for n in range(n):

        last_moments = []

        for r in range(curve.fractal):
            current_moment = []
    
            for k in range(curve.genus):
        
                new_moment = moment_sub_k[int(curve.base_maps[r][k][0])]
        
                if curve.base_maps[r][k][-1]=='1':
                    new_moment = np.flip(den*curve.genus**n - new_moment)
        
                norm_moment = new_moment + k*den*curve.genus**n
                current_moment = np.concatenate((current_moment, norm_moment))
            
            last_moments.append(current_moment)
    
        moment_sub_k = last_moments
        
    moment_sub_k.append(den)
        
    return [int32(k) for k in moment_sub_k]


def get_junction_ratio(junction, n):
       
    bm1,bm2 = junction
    
    # Генерируем моменты n-ого подразделения
    moments = curve.get_vertex_moments()
    moment_sub_k = get_moment_sub_k(moments,n)
    
    first_moments = moment_sub_k[int(bm1[0])]
    second_moments = moment_sub_k[int(bm2[0])]
    
    if bm1[-1] == '1':
        first_moments = np.flip(moment_sub_k[-1]*curve.genus**n - first_moments)
        
    if bm2[-1] == '1':
        second_moments = np.flip(moment_sub_k[-1]*curve.genus**n - second_moments)
    
    d = np.arange(0, 2**curve.dim*curve.genus**(n), 2**curve.dim)
    d = np.flip(d)
        
    first_moments = np.delete(first_moments, d, 0)
    second_moments = np.delete(second_moments, d, 0)
    
    if n != 0:    
        second_moments = moment_sub_k[-1]*curve.genus**n + second_moments
    else:
        second_moments = moment_sub_k[-1] + second_moments

    # Генерируем n-ое подразделение кривой
    sub = curve.get_subdiv(n,plot=False)

    frac_1 = curve.get_fraction(sub[int(bm1[0])],bm1[1:])
    coord_1 = curve.get_curve_coord(frac_1)

    frac_2 = curve.get_fraction(sub[int(bm2[0])],bm2[1:])
    coord_2 = curve.get_curve_coord(frac_2,start=coord_1[-1])

    coord_1.pop(0)
    coord_2.pop(0)
    
    coord_1 = np.array(coord_1)
    coord_2 = np.array(coord_2)

    x1 = coord_1
    x2 = coord_2

    t1 = first_moments
    t2 = second_moments
    
    if curve.dim == 2:

        r = get_ratio_l2_2d(x1,x2,t1,t2)
        
    elif curve.dim == 3:
        
        r = get_ratio_l2_3d(x1,x2,t1,t2)
        
    return moment_sub_k[-1]*r




junctions = curve.get_junctions()

t = time.time()
ratio = [get_junction_ratio(k, sub_numb) for k in junctions]
elapsed = time.time() - t
print('время нахождения всех отношений -',elapsed)

for k in range(len(ratio)):
    print(junctions[k],'-',ratio[k])

