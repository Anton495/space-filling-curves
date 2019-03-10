from examples import *
from pl_curves import plcurves

#Indicate subdivision number,curve number and run
subdiv_number = 2
curve_number = 1

get_Curve =[#2d curve examples (mono-fractal)
            get_Hilbert_Curve(),         # 0
            get_Peano_Curve(),           # 1
            get_Meurthe_Curve_Curve(),   # 2
            get_Coil_Curve_Curve(),      # 3
            get_Serpentine_Curve(),      # 4
            get_R_Curve(),               # 5
            #2d curve examples (bi-fractal)
            get_Beta_Omega_Curve(),      # 6
            #2d curve examples (tetra-fractal)
            get_ARW_Curve(),             # 7
            #3d curve examples (mono-fractal)
            get_Tokarev_Curve(),         # 8 
            get_Haverkort_Curve_1(),     # 9
            get_Haverkort_Curve_2(),     # 10
            #3d curve examples (bi-fractal)
            get_Neptunus_Curve(),        # 11
            get_Luna_Curve()]            # 12

curve = get_Curve[curve_number]

plcurves(curve.get_subdiv_n(subdiv_number), curve.dim, curve.genus, subdiv_number)
