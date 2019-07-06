from examples import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_curve(sub_n, dim, genus, n):
    
    # Масштабируем координаты, делим все значения на genus**((n+1)/dim)
    sub_n = [[i/(genus**((n+1)/dim)) for i in sub_n[j]] for j in range(len(sub_n))]   
    # Находим минимум по каждой коодинате x,y,z b и т.д.,
    # это необходимо для кривых, у которых точки входа и выхода находятся не в вершинах 
    min_col = [min([row[i] for row in sub_n]) for i in range(dim)]
    # Сдвиг всех координат на 1/(2*genus**((n+1)/(dim))) с учетом минимальных значений по координатам
    sub_n = [[sub_n[j][i] - min_col[i] + 1/(2*genus**((n+1)/(dim))) 
                 for i in range(dim)] for j in range(len(sub_n))]
    # Cетка для графика
    ticks = list(range(0,int(genus**(n/2)+1)))
    ticks = [i/ticks[-1] for i in ticks]
    
    if dim == 2:
       
        fig = plt.figure()
        ax = fig.gca()
        ax.tick_params(direction='in', top=True, right=True)
        plt.gcf().set_size_inches(9,9)        
        plt.xticks(ticks,[])
        plt.yticks(ticks,[])
        plt.grid()
        plt.axis([0,1,0,1])
        sub1 = [row[0] for row in sub_n]
        sub2 = [row[1] for row in sub_n]
        plt.plot(sub1,sub2,'k')

    elif dim == 3:
    
        fig = plt.figure(figsize = (9,9))
        ax = fig.gca(projection='3d')
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        ax.set_zlim(0,1)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_zticks(ticks)
        ax.tick_params(colors = 'w')
        sub1 = [row[0] for row in sub_n]
        sub2 = [row[1] for row in sub_n]
        sub3 = [row[2] for row in sub_n]
        ax.plot(sub1,sub2,sub3,'k')
        
    plt.savefig('curve_image.png')
