from examples import *
from pl_curves import plcurves

#Indicate subdivision number and run
k = 3

get_Curve =[#2d curve examples (mono-fractal)
            get_Hilbert_Curve(k),        # 0
            get_Peano_Curve(k),          # 1
            get_Meurthe_Curve(k),  # 2
            get_Coil_Curve(k),     # 3
            get_Serpentine_Curve(k),     # 4
            get_R_Curve(k),              # 5
            #2d curve examples (bi-fractal)
            get_Beta_Omega_Curve(k),     # 6
            #2d curve examples (tetra-fractal)
            get_ARW_Curve(k),            # 7
            #3d curve examples (mono-fractal)
            get_Tokarev_Curve(k),        # 8 
            get_Haverkort_Curve_1(k),    # 9
            get_Haverkort_Curve_2(k),    # 10
            #3d curve examples (bi-fractal)
            get_Neptunus_Curve(k),       # 11
            get_Luna_Curve(k)]           # 12

#Indicate the curve number and run
curve = get_Curve[0]

plcurves(curve.subdiv_n, curve.dim, curve.genus, curve.number_sub)
