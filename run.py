from examples import *
from pl_curves import plot_curve

#Indicate subdivision number,curve number and run
subdiv_number = 4
curve_number = 7

get_curve =[#2d curve examples (mono-fractal)
            get_hilbert_curve(),         # 0
            get_peano_curve(),           # 1
            get_meurthe_curve(),         # 2
            get_coil_curve(),            # 3
            get_serpentine_curve(),      # 4
            get_R_curve(),               # 5
            #2d curve examples (bi-fractal)
            get_beta_Omega_curve(),      # 6
            #2d curve examples (tetra-fractal)
            get_ARW_curve(),             # 7
            #3d curve examples (mono-fractal)
            get_tokarev_curve(),         # 8 
            get_haverkort_curve_1(),     # 9
            get_haverkort_curve_2(),     # 10
            #3d curve examples (bi-fractal)
            get_neptunus_curve(),        # 11
            get_luna_curve()]            # 12

curve = get_curve[curve_number]

plot_curve(curve.get_subdiv(subdiv_number), curve.dim, curve.genus, subdiv_number)
