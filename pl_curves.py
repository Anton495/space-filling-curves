import matplotlib.pyplot as plt
from numpy import linspace
from mpl_toolkits.mplot3d import Axes3D

def plcurves(subdiv_n, dim, genus, number_sub):

    if dim == 2:

        ticks = linspace(0,1,genus**(number_sub/2)+1)
    
        plt.figure()
        plt.gcf().set_size_inches(9,9)
        plt.xticks(ticks,[])
        plt.yticks(ticks,[])
        plt.grid(True)
        plt.axis([0,1,0,1])
        plt.plot(subdiv_n[:,0],subdiv_n[:,1],'k')

    elif dim == 3:
    
        ticks = linspace(0,1,genus**(number_sub/3)+1)
    
        fig = plt.figure(figsize = (9,9))
        ax = fig.gca(projection='3d')
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        ax.set_zlim(0,1)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_zticks(ticks)
        ax.tick_params(colors = 'w')
        ax.plot(subdiv_n[:,0],subdiv_n[:,1],subdiv_n[:,2],'k')
    

    
    
    
    
    
    
    