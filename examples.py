from space_filling_curve import FractalCurve

#2d curve examples (mono-fractal)

def get_Hilbert_Curve(k):
    return FractalCurve(
            number_sub = k,
            alphabet = 'iIjJ',
            proto = [['jiJ']],
            base_maps = [['ji','ij','ij','JI']]
    )
    
def get_Peano_Curve(k):
    return FractalCurve(
            number_sub = k,
            alphabet = 'iIjJ',
            proto = [['jjiJJijj']],
            base_maps = [['ij','Ij','ij','iJ','IJ','iJ','ij','Ij','ij']]
    )
    
def get_Meurthe_Curve_Curve(k):
    return FractalCurve(
            number_sub = k,
            alphabet = 'iIjJ',
            proto = [['jjiJJijj']],
            base_maps = [['ij','Ji','ij','jI','JI','iJ','ji','Ji','ij']]
    )
    
def get_Coil_Curve_Curve(k):
    return FractalCurve(
            number_sub = k,
            alphabet = 'iIjJ',
            proto = [['jjiJJijj']],
            base_maps = [['ji','Ji','ji','jI','JI','jI','ji','Ji','ji']]
    ) 
    
def get_Serpentine_Curve(k):
    return FractalCurve(
            number_sub = k,
            alphabet = 'iIjJ',
            proto = [['jjiJJijj']],
            base_maps = [['ij','Ji','ji','iJ','JI','iJ','ji','Ji','ij']]
    )    
    
def get_R_Curve(k):
    return FractalCurve(
            number_sub = k,
            alphabet = 'iIjJ',
            proto = [['jjiiJIJi']],
            base_maps = [['ji','ji','ij','ij','ij','IJ','JI','JI','ij']]
    )    
    
#2d curve examples (bi-fractal)
    
def get_Beta_Omega_Curve(k):
    return FractalCurve(
            number_sub = k,
            alphabet = 'iIjJ',
            proto = [['jiJ'],
                     ['jiJ']],
            base_maps = [['1iJ0','1Ji0','1ji1','1IJ1'],
                         ['1iJ0','1Ji0','1ji1','0jI0']]
    )    
    
#2d curve examples (quater-fractal)
    
def get_ARW_Curve(k):
    return FractalCurve(
            number_sub = k,
            alphabet = 'iIjJ',
            proto = [['i','Ij','i'],
                     ['jiJ'],
                     ['jiJ'],
                     ['jiJ']],
            base_maps = [['3ij0','1jI1','2Ji0','1iJ0'],
                         ['3ji0','2Ij1','1ij0','1JI0'],
                         ['0ji0','1Ji0','0jI0','1JI0'],
                         ['0ij0','2Ji0','0jI0','3Ji1']]
    ) 
    
#3d curve examples (mono-fractal)

def get_Tokarev_Curve(k):
    return FractalCurve(
            number_sub = k,
            alphabet = 'iIjJkK',
            proto = [['kjKikJK']],
            base_maps = [['jki','kij','kij','iJK','iJK','KIj','KIj','JkI']]
    ) 

def get_Haverkort_Curve_1(k):
    return FractalCurve(
            number_sub = k,
            alphabet = 'iIjJkK',
            proto = [['kjKikJK']],
            base_maps = [['kji0','jik0','kIj1','iKJ0','IKJ1','KIj0','Kij1','Jki1']]
    ) 

def get_Haverkort_Curve_2(k):
    return FractalCurve(
            number_sub = k,
            alphabet = 'iIjJkK',
            proto = [['kjKikJK']],
            base_maps = [['KIJ1','KJI1','KjI0','Jki1','jki0','kjI1','kJI0','iKJ0']]
    ) 

#3d curve examples (bi-fractal, Haverkort curves)
    
def get_Neptunus_Curve(k):
    return FractalCurve(
            number_sub = k,
            alphabet = 'iIjJkK',
            proto = [['kjKikJK'],
                     ['kiKjIki']],
            base_maps = [['0kji','1kji','1KiJ','1jKI','1ikj','1KJi','0kJI','1jKI'],
                         ['0jki','1jki','1iKJ','0KiJ','1JiK','1IKj','0ikj','1ijk']]
    )

def get_Luna_Curve(k):
    return FractalCurve(
            number_sub = k,
            alphabet = 'iIjJkK',
            proto = [['kjKikJK'],
                     ['kiKjIki']],
            base_maps = [['1ijk','0KJi','1KiJ','1jKI','1jik','1IKj','0kJI','1kJI'],
                         ['1jik','0JKi','1iKJ','0KiJ','1KjI','1JIk','0ikj','1ikj']]
    )    
