import pandas as pd

from atmospere_utils import *

###################################
# classes / utils for monte carlo #
###################################

class RandomBlackBody:
    def __init__(self, T, _n, fname):
        self.T = T
        self.n = _n

        f_table = pd.read_csv(fname, names=['nlT', 'f'])
        self.nlT = f_table['nlT'].to_numpy()
        self.f = f_table['f'].to_numpy()

    def get_lambda(self) -> float:
        """get random wavelength for black surface
        
        return: lambda - float [mum]
        """
        R_lambda = np.random.random()

        # linearly interpolate black body emissive power table
        return f_lerp(R_lambda, self.f, self.nlT) / (self.T * self.n)

random_rectangle = lambda W, L: np.array([(np.random.random() - .5) * W, (np.random.random() - .5) * L])
"""get random location on rectange

return: (x, y) - (float, float) [m]
"""

random_range = lambda start, stop: (stop - start) * np.random.random() + start
"""get random location on range

return: value - float
"""

random_phi = lambda: 2 * np.pi * np.random.random()
"""get random azimuthal angle with azimuthal symmetry

return: phi - float [rad]
"""

theta_rayleigh = np.linspace(0.0, np.pi, num=150)
R_theta_rayleigh = (1 - .75 * (np.cos(theta_rayleigh)**3 / 3 + np.cos(theta_rayleigh))) / 2
random_theta_rayleigh = lambda: f_lerp(np.random.random(), R_theta_rayleigh, theta_rayleigh)
"""get random polar angle for rayleigh scattering

return: theta - float [rad]
"""

random_theta_emit = lambda: np.arccos(2 * np.random.random() - 1)
"""get random polar angle for reemmission after absorption

return: theta - float [rad]
"""

random_l_s = lambda _lambda: np.log(1 / np.random.random()) / get_sigma_s(_lambda)
"""get random scattering length

return: l_scattering - float [m]
"""

random_l_a = lambda _lambda: np.log(1 / np.random.random()) / get_kappa(_lambda)
"""get random absorption length

return: l_absorption - float [m]
"""