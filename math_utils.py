import numpy as np
from datetime import datetime
from constants import L_AIR

##############
# math utils #
##############

# angle conversions
r2d = 180. / np.pi
d2r = np.pi / 180.

quad_sqrt = lambda d: np.sqrt(d) if d >= 0 else np.nan

# quadratic formula
quad = lambda a, b, c: ((-b + quad_sqrt(b**2 - 4 * a * c)) / (2 * a), (-b - quad_sqrt(b**2 - 4 * a * c)) / (2 * a))
"""everyone's favorite quadratic formula  

    args:  
        a - float: quadratic coefficient   
        b - float: linear coefficient   
        c - float: constant coefficient   

    returns: (plus, minus) - tuple[float, float]: plus sqrt and minus sqrt
"""

# linear interpolation
lerp = lambda x, x_0, y_0, x_1, y_1: y_0 + (x - x_0) * (y_1 - y_0) / (x_1 - x_0)

counter_clockwise = lambda theta: np.array([
    [np.cos(theta), -np.sin(theta)], 
    [np.sin(theta),  np.cos(theta)]
])

clockwise = lambda theta: np.array([
    [ np.cos(theta), np.sin(theta)], 
    [-np.sin(theta), np.cos(theta)]
])

xify        = lambda arr: np.array([   0.0, *arr, np.inf])
functionify = lambda arr: np.array([arr[0], *arr, arr[-1]])

def f_lerp(_x, x: np.array, f: np.array):
    """linear interpolate f(x) if map x -> f is expressed as 2 sorted np.arrays

    args:
        _x: x value used in interpolation of f(x)
            ie f(_x) therefore same datatype and within the bounds of x
        x - np.array: sorted x values
        f - np.array: f(x_i) -> f_i
    """
    big_idx = np.searchsorted(x, _x)
    if big_idx == len(x): return f[-1]

    small_idx = big_idx - 1

    # linearly interpolate
    return lerp(_x, x[small_idx], f[small_idx], x[big_idx], f[big_idx])

# print(f'ccw : {counter_clockwise(np.pi/6)}')
# print(f'cw : {clockwise(np.pi/6)}')

sph2rec = lambda r, phi, theta: np.array([
    r * np.cos(phi) * np.sin(theta),
    r * np.sin(phi) * np.sin(theta),
    r * np.cos(theta)
])
"""spherical to rectangular coordinates  
args:  
    r, phi, theta - cylindrical coordinates, radius, azimuthal and polar angles  
returns:  
    x, y, z - rectangular coordinates  
"""

rec2sph = lambda x, y , z: np.array([
    r := np.sqrt(x**2 + y**2 + z**2),
    np.atan2(y, x),
    np.atan2(z, r)
])
"""rectangular to spherical coordinates  
args:  
    x, y, z - rectangular coordinates  
returns:  
    r, phi, theta - cylindrical coordinates, radius, azimuthal and polar angles  
"""


def rayleigh_airmass(altitude: float, zenith_length: float = L_AIR) -> float:
    """caclulate the length of air dense enough for relevent rayleigh scattering  

    [3] K. Pickering, “Reidentifcation of some entries in the Ancient Star 
        Catalog” The International Journal of Scientific History. vol. 12, 
        Sep 2002.  

    h is apparent altitude in degrees  
    1 / sin(h + 244 / (165 + 47 * h^1.1))  

    args:  
        altitude - float: apparent altitude in radians  
        zenith_length - float: rayleigh scattering length at zenith, default L_AIR = 50 [km]  
    """
    alt_d = r2d * altitude
    h_rad = d2r * (alt_d + 244 / (165 + 47 * alt_d**1.1))
    return zenith_length / np.sin(h_rad)

fmt_now = lambda: datetime.now().strftime('%H:%M:%S')