import numpy as np

######################
# physical constants #
######################

c_0   = 299792458        # speed of light in vacuum [m/s]
h     = 6.62607015e-34   # Planck constant [J s]
k     = 1.380649e-23     # Boltzmann constant [J/K]
q     = 1.6021766208e-19 # electron charge [C]
sigma = 5.670374419e-8   # Stefan-Boltzmann constant [W/(m^2 K^4)]
n_air = 1                # assumed refractive index of air

########################
# simulation constants #
########################

L_AIR    = 50_000 # [m]
W_T, L_T = 30, 200 # [m]
GROUND   = -15 # [m]
M_PARABOLA = np.array([.06, .033, .025, .0175])
RHO_PARABOLA = np.array([1.0, .98, .96, .94])
ALPHA_TUBE = np.array([1.0, .98, .96, .94])
N_RAYS = 10_000