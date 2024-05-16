# originally made with more ambitious ideas in mind
# all we need is the base_fct for our code

import numpy as np

class _2d_wavelet:
    def __init__(self, name):
        self.name = name
        if name == 'mexh':
            self.cpsi = 2. * np.pi**2.
        
    def base_fct(self, x):
        r_sq = np.linalg.norm(x, axis = -1)**2.
        return ( 2 - r_sq ) * np.exp(-r_sq / 2.)
    
    def fct(self, x, b, a):
        arg = (x - b) / a
        r_sq = np.linalg.norm(arg, axis = 1)**2.
        return ( 2 - r_sq ) * np.exp(-r_sq / 2.)
    
    def estimate_coefficients(self, x, b, a):
        fct_array = self.fct(x, b, a)
        return np.sum(fct_array) / a / len(x)

