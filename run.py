from examples import *
from pl_curves import plcurves

#Indicate subdivision number 
k = 3

#Choose the curve

#2d curve examples (mono-fractal)
curve = get_Hilbert_Curve(k)
#curve = get_Peano_Curve(k)
#curve = get_Meurthe_Curve_Curve(k)
#curve = get_Coil_Curve_Curve(k)
#curve = get_Serpentine_Curve(k)
#curve = get_R_Curve(k)

#2d curve examples (bi-fractal)
#curve = get_Beta_Omega_Curve(k)

#2d curve examples (quater-fractal)
#curve = get_ARW_Curve(k)

#3d curve examples (mono-fractal)
#curve = get_Tokarev_Curve(k)
#curve = get_Haverkort_Curve_1(k)
#curve = get_Haverkort_Curve_2(k)

#3d curve examples (bi-fractal)
#curve = get_Neptunus_Curve(k)
#curve = get_Luna_Curve(k)


plcurves(curve.subdiv_n, curve.dim, curve.genus, curve.number_sub)
