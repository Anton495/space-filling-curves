from space_filling_curve import FractalCurve
from convert_bm import convert_bm

#2d curve examples (mono-fractal)

def get_hilbert_curve():
    return FractalCurve(
            chain_proto = ['jiJ'],
            base_maps = [['ji','ij','ij','JI']]
    )
    
def get_peano_curve():
    return FractalCurve(
            chain_proto = ['jjiJJijj'],
            base_maps = [['ij','Ij','ij','iJ','IJ','iJ','ij','Ij','ij']]
    )
    
def get_meurthe_curve():
    return FractalCurve(
            chain_proto = ['jjiJJijj'],
            base_maps = [['ji','Ji','ij','jI','JI','iJ','ji','Ji','ij']]
    )

def get_coil_curve():
    return FractalCurve(
            chain_proto = ['jjiJJijj'],
            base_maps = [['ji','Ji','ji','jI','JI','jI','ji','Ji','ji']]
    )
    
def get_serpentine_curve():
    return FractalCurve(
            chain_proto = ['jjiJJijj'],
            base_maps = [['ij','Ji','ji','iJ','JI','iJ','ji','Ji','ij']]
    )    
    
def get_R_curve():
    return FractalCurve(
            chain_proto = ['jjiiJIJi'],
            base_maps = [['ji','ji','ij','ij','ij','IJ','JI','JI','ij']]
    )    
    
#2d curve examples (bi-fractal)
    
def get_beta_Omega_curve():
    return FractalCurve(
            chain_proto = ['jiJ','jiJ'],
            base_maps = [['1iJ','1Ji','1ji1','1IJ1'],
                         ['1iJ','1Ji','1ji1','0jI' ]]
    )    
    
#2d curve examples (quater-fractal)
    
def get_ARW_curve():
    return FractalCurve(
            chain_proto = [['i','Ij','i'],'jiJ','jiJ','jiJ'],
            base_maps = [['3ij','1jI1','2Ji','1iJ' ],
                         ['3ji','2Ij1','1ij','1JI' ],
                         ['0ji','1Ji' ,'0jI','1JI' ],
                         ['0ij','2Ji' ,'0jI','3Ji1']]
    )
    
#3d curve examples (mono-fractal)

def get_tokarev_curve():
    return FractalCurve(
            chain_proto = ['kjKikJK'],
            base_maps = [['jki','kij','kij','iJK','iJK','KIj','KIj','JkI']]
    ) 

def get_ye_curve():
    return FractalCurve(
            chain_proto = ['jiJijjIIjjiJijiiJIJiJIJi'],
            base_maps = [['jI1','ij','ij','Ji1','jI1',
                          'jI1','jI1','IJ','IJ','jI1',
                          'ij','ij','Ji1','jI1','ij',
                          'ij','ij','IJ','Ji1','Ji1',
                          'ij','IJ','Ji1','Ji1','ij']]
    ) 
    
def get_peano7_curve():
    return FractalCurve(
            chain_proto = ['jiJijiJijjIjijIIJJIIjijIjjiJijiJijiiJIJiJIJiJIJi'],
            base_maps = [['jI1','ij','ij','Ji1','jI1','ij','ij',
                           'Ji1','jI1','jI1','jI1','IJ','ij','jI1',
                           'jI1','IJ','Ji1','Ji1','IJ','IJ','IJ',
                           'ij','jI1','jI1','IJ','jI1','ij','ij',
                           'Ji1','jI1','ij','ij','Ji1','jI1','ij',
                           'ij','ij','IJ','Ji1','Ji1','ij','IJ',
                           'Ji1','Ji1','ij','IJ','Ji1','Ji1','ij']]
    )
    
def get_bauman_curve():
    return FractalCurve(
            chain_proto = ['ijIjiiJJ'],
            base_maps = [['ij','ji','jI1','iJ1','ij','ij','Ij1','Ji1','Ji1']]
    )

def get_haverkort_curve_1():
    return FractalCurve(
            chain_proto = ['kjKikJK'],
            base_maps = [['kji','kij','kIj1','iKJ','IKJ1','KIj','Kij1','Jki1']]
    ) 

def get_haverkort_curve_2():
    return FractalCurve(
            chain_proto = ['kjKikJK'],
            base_maps = [['KIJ1','KJI1','KjI','Jki1','jki','kjI1','kJI','iKJ']]
    ) 

#3d curve examples (bi-fractal, Haverkort curves)
    
def get_neptunus_curve():
    return FractalCurve(
            chain_proto = ['kjKikJK','kiKjIki'],
            base_maps = [['0kji','1kji','1KiJ','1jKI','1ikj','1KJi','0kJI','1jKI'],
                         ['0jki','1jki','1iKJ','0KiJ','1JiK','1IKj','0ikj','1ijk']]
    )

def get_luna_curve():
    return FractalCurve(
            chain_proto = ['kjKikJK','kiKjIki'],
            base_maps = [['1ijk','0KJi','1KiJ','1jKI','1jik','1IKj','0kJI','1kJI'],
                         ['1jik','0JKi','1iKJ','0KiJ','1KjI','1JIk','0ikj','1ikj']]
    )    

#3d curve examples (bi-fractal, new curves)

def get_17_curve():
    return FractalCurve(
            chain_proto = ['jkiKJkI','jiJkjIJ'],
            base_maps = convert_bm([['1JKI~','0jKI','1kji','0kiJ~','1KiJ','0JKi','1kJI','1IkJ'],
                                    ['1JKI~','0Ijk~','0jiK~','0KJI~','0Jki~','0ijK~','0IjK','1JIk']])

