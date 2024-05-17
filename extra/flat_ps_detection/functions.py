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
    
def closest_grid_point_(point, grid_points):
    relative_positions = grid_points - point
    distances = np.linalg.norm(relative_positions, axis = -1)
    bxby_index_shape = distances.shape
    flattened_index = np.argmin(distances)
    return np.unravel_index(flattened_index, bxby_index_shape)

