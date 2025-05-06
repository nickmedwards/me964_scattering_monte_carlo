from constants import *
from math_utils import *

from mc_utils import RandomBlackBody, random_range
import numpy.linalg as linalg

########################
# sky and ray geometry #
########################

class Ray:
    def __init__(self, lam: float, origin: np.array, hat: np.array):
        self.lam = lam
        self.origin = origin
        self.hat = hat
        
    def __repr__(self):
        return f'\n\nlambda:\t{self.lam}\norigin:\t{self.origin}\nhat:\t{self.hat}\n'
    
    def __call__(self, t) -> np.array:
        return np.array(self.origin + t * self.hat)

get_l_ground = lambda ray, gr = GROUND: (gr - ray.origin[2]) / ray.hat[2]

class Sky:
    """an approximation of the amount of rayleigh scattering air based on z position"""
    N_POINTS = 181

    def __init__(self, L_air: float = L_AIR, ground: float = GROUND):
        """
        args:  
            L_air - float: modeled length of air at zenith, default L_AIR = 50 [km]
            ground - float: small offset to allow tube to have center at and troughs to have foci at x,z -> 0,0, default = -15 [m]
        """
        self.z_max = L_air
        self.z_min = ground
        self.z = np.linspace(ground, L_air, num=Sky.N_POINTS)
        self.altitude = np.linspace(0.0, np.pi / 2, num=Sky.N_POINTS)
        self.air = np.array([rayleigh_airmass(a, L_air) for a in self.altitude])
        self.flipped_air = np.flip(self.air) # flip for searching
    
    def get_l_sky(self, ray: Ray):
        # approximate height with xz and xy plane intercepts
        z = ray.origin[2]
        # if higher than atmosphere then force absorption in hit_trace
        if z > self.z_max: return -np.inf

        # approximate global altitude angle with ray polar angle
        approx_alt = np.pi / 2 - rec2sph(*ray.hat)[2]

        # print(f'z_step: {z_step}')
        # print(f'possible_idx: {possible_idx}')
        # print(f'small_idx: {small_idx}')
        # print(f'big_idx: {big_idx}')
        # print(f'self.air[small_idx]: {self.air[small_idx]}')
        # print(f'self.air[big_idx]: {self.air[big_idx]}')
        # print(f'plane_airmass: {plane_airmass}')

        # airmass at that z intercept
        # plane_airmass = f_lerp(z, self.z, self.air)
        # airmass_idx = len(self.air) - np.searchsorted(self.flipped_air, plane_airmass)
        # since z is linearly spaced can skip f_lerp
        # z_step = (self.z_max - self.z_min) / (Sky.N_POINTS - 1)
        possible_z_idx = (Sky.N_POINTS - 1) * (z - self.z_min) / (self.z_max - self.z_min)
        airmass_idx = int(np.ceil(possible_z_idx))

        temp_points = Sky.N_POINTS - airmass_idx - 1
        # print(f'z: {z}')
        # print(f'possible_z_idx: {possible_z_idx}')
        # print(f'airmass_idx: {airmass_idx}')
        # print(f'temp_points: {temp_points}')

        if temp_points == 0: return self.air[-1] - z
        temp_airmass = self.air[airmass_idx:]

        # interpolate rescaled airmass function
        alt_step = (np.pi / 2) / (Sky.N_POINTS - airmass_idx - 1)
        possible_idx = approx_alt / alt_step
        small_idx, big_idx = int(np.floor(possible_idx)), int(np.ceil(possible_idx))
        
        # print(f'approx_alt: {approx_alt}')
        # print(f'airmass_idx: {airmass_idx}')
        
        # print(f'alt_step: {alt_step}')
        # print(f'len: {(Sky.N_POINTS - airmass_idx - 1)}')
        # print(f'alt_step * len: {alt_step * (Sky.N_POINTS - airmass_idx - 1)}')
        # print(f'possible_idx: {possible_idx}')
        # print(f'small_idx: {small_idx}')
        # print(f'big_idx: {big_idx}')
        approx_airmass = temp_airmass[small_idx] if small_idx == big_idx else lerp(approx_alt, small_idx * alt_step, temp_airmass[small_idx], big_idx * alt_step, temp_airmass[big_idx])
        # print(f'approx_airmass: {approx_airmass}')
        # print(f'return? {approx_airmass - z}')
        
        return approx_airmass - z
        # # angles[1] = (np.pi / 2) / (Sky.N_POINTS - airmass_idx - 1)
        # # (np.pi / 3.1) / angles[1] -> posible index
        # # if floor == ciel: possibe is index, otherwise lerp(angle, floor * angles[1], self.air[airmass_idx:][floor], ceil * angles[1], self.air[airmass_idx:][ceil])
        # return f_lerp(
        #     np.pi / 2 - hat_sph[2], 
        #     np.array([*np.linspace(0, np.pi / 2, num=Sky.N_POINTS - airmass_idx), np.inf]),
        #     np.array([*self.air[airmass_idx:], self.air[-1]])
        # ) - z

#####################
# hittable geometry #
#####################

def reflect_2d(n_hat: np.array, _r_hat: np.array) -> np.array:
    # rescale ray.hat for xz plane
    r_xz = np.array([_r_hat[0], _r_hat[2]])
    r_hat = r_xz / linalg.norm(r_xz)
    # (-r_hat) dot n_hat = cos theta
    # theta = acos(-r_hat_x * n_hat_x - r_hat_z * n_hat_z)
    # sometimes bc of float math, 
    # unit vectors w unit norm have above unit args to acos
    # r:   [-0.70746101  0.70675238]  1.0
    # n:   [ 0.70746101 -0.70675238]  1.0
    # arg: 1.0000000000000002
    # so buffer
    arg = -r_hat[0] * n_hat[0] - r_hat[1] * n_hat[1]
    arg =  1 if arg >  1 else arg
    arg = -1 if arg < -1 else arg
    theta_xz = np.acos(arg)

    # rotates ccw if ray has positive z component in n_hat coords 
    n_theta_xz = np.atan2(n_hat[1], n_hat[0])
    rotate = counter_clockwise if (clockwise(n_theta_xz) @ np.array([r_hat[0], r_hat[1]]))[1] > 0 else clockwise

    # print(f'n_hat:\t{n_hat}')
    # print(f'theta_xz:\t{r2d*theta_xz}')
    # print(f'n_theta_xz:\t{r2d*n_theta_xz}')
    # print(f'test = {clockwise(n_theta_xz) @ np.array([r_hat[0], r_hat[1]])}')
    # print('counter_clockwise' if (clockwise(n_theta_xz) @ np.array([r_hat[0], r_hat[1]]))[1] > 0 else 'clockwise')
    # new_r_hat_xz = -rotate(2 * theta_xz) @ np.array([r_hat[0], r_hat[1]])
    # print(f'new_r_hat_xz:\t{new_r_hat_xz}')

    return -rotate(2 * theta_xz) @ np.array([r_hat[0], r_hat[1]])

class Trough:
    N_POINTS = 31
    QUADRATIC_THRESHOLD = 3 * 10**-9

    def __init__(self, W: float, L: float, m: float, rho: float):
        """parabolic trough geometry with focus at origin of x, z plane  
        parabolic parameter, m forces focus to be at origin  
        surface described ->  
            -W/2 <= x <= W/2  
            -L/2 <= y <= L/2  
            z = m * x^2 - 1 / (4 * m)  

        args:  
            W   - float: trough width  
            L   - float: trough length  
            m   - float: parabolic parameter  
            rho - float: reflectance  
        """
        self.W = W
        self.L = L
        self.m = m
        self.rho = rho
        # s = Integrate[Sqrt[1 + (2 m x)^2], x] evaulated from -W/2 to W/2
        # s = 1/2 x sqrt(4 m^2 x^2 + 1) - log(sqrt(4 m^2 x^2 + 1) - 2 sqrt(m^2) x)/(4 sqrt(m^2)) + constant
        # symmetric so evalulate from 0 to W/2 and multiply by 2
        # 2*m*x at x = W/2 -> m*W
        half_s = (W/2) * np.sqrt((m * W)**2 + 1) / 2 - np.log(np.sqrt((m * W)**2 + 1) - m * W) / (2 * m)
        self.surface_area = L * 2 * half_s

        self.z_t = m * (W/2)**2 - 1 / (4 * m)
        self.z_min = - 1 / (4 * m)
        self.H_t = self.z_t - self.z_min

    def for_plot(self):
        x = np.linspace(-self.W/2, self.W/2, num=Trough.N_POINTS)
        y = np.linspace(-self.L/2, self.L/2, num=Trough.N_POINTS)
        z = self.m * x**2 - 1 / (4 * self.m)
        plot_x, plot_y = np.meshgrid(x, y)
        return plot_x, plot_y, z
    
    def enters(self, x, y):
        return (np.abs(x) <= self.W / 2) and (np.abs(y) <= self.L / 2)
    
    def hit_bounds(self, x, y, z):
        return (np.abs(x) <= self.W / 2) and (np.abs(y) <= self.L / 2) and (z < self.z_t) and (z > self.z_min)
    
    def is_hit(self, ray: Ray):
        # if ray is barely moving in x
        # simplify solution bc reduces to linear equation
        if ray.hat[0] < Trough.QUADRATIC_THRESHOLD:
            t = (self.m * ray.origin[0]**2 - 1 / (4 * self.m) - ray.origin[2]) / ray.hat[2]
            return (self.hit_bounds(*(hit_point := ray(t))), hit_point, t) if t > 0 else (False, np.array([np.nan, np.nan, np.nan]), 0)
        
        a = self.m * ray.hat[0]**2
        b = (2 * self.m * ray.origin[0] * ray.hat[0] - ray.hat[2])
        c = self.m * ray.origin[0]**2 - 1 / (4 * self.m) - ray.origin[2]
        # enforce t > 0
        # parabola only hittable if split (t_m, t_p) -> (+, -) or (-, +)
        t_p, t_m = quad(a, b, c)
        # print(f'a: {a}\tb: {b}\tc: {c}\t')
        # print(f't_m: {t_m}\tt_p: {t_p}')
        if t_p * t_m < 0:
            t = t_p if t_p > 0 else t_m
            hit_point = ray(t)
            return (self.hit_bounds(*hit_point), hit_point, t)
        return (False, np.array([np.nan, np.nan, np.nan]), 0)

    def reflect(self, ray: Ray, point: np.array):
        # print('TROUGH REFLECT')
        # dz/dx = 2 * m * x -> (1 , 2 * m * x)
        # normal to dz/dx -> (-2 * m * x, 1)
        n_dz_dx = np.array([-2 * self.m * point[0], 1])
        n_hat = n_dz_dx / linalg.norm(n_dz_dx)
        new_r_hat_xz = reflect_2d(n_hat, ray.hat)
        # dz_dx = np.array([1, 2 * self.m * point[0]])
        # print(f'dz_dx:{dz_dx}')
        # d_hat = dz_dx / linalg.norm(dz_dx)
        # print(f'd_hat:{d_hat}')
        # print(f'n_hat:{n_hat}')
        # print(f'new_r_hat_xz:{new_r_hat_xz}')
        new = np.array([new_r_hat_xz[0], ray.hat[1], new_r_hat_xz[1]])
        r_hat = new / linalg.norm(new)
        return Ray(ray.lam, point, r_hat)

class Tube:
    N_POINTS = 30
    HIT_THRESHOLD = 10**-9

    def __init__(self, L: float, r: float, alpha: float):
        """tube centered at origin of x, z plane  

        surface described ->  
            x^2 + z^2 = r^2  
            -L/2 <= y <= L/2  

        args:  
            L     - float: tube length  
            r     - float: tube radius  
            alpha - float: absorptance  
        """
        self.L = L
        self.r = r
        self.alpha = alpha
        self.surface_area = 2 * np.pi * r * L

    def for_plot(self):
        x = np.array([self.r * np.cos(theta) for theta in np.linspace(0, 2 * np.pi, num=Tube.N_POINTS)])
        y = np.linspace(-self.L/2, self.L/2, num=Tube.N_POINTS)
        z = np.array([self.r * np.sin(theta) for theta in np.linspace(0, 2 * np.pi, num=Tube.N_POINTS)])
        plot_x, plot_y = np.meshgrid(x, y)
        return plot_x, plot_y, z
    
    def is_hit(self, ray: Ray):
        a = ray.hat[0]**2 + ray.hat[2]**2
        b = 2 * (ray.origin[0] * ray.hat[0] + ray.origin[2] * ray.hat[2])
        c = ray.origin[0]**2 + ray.origin[2]**2 - self.r**2

        t_p, t_m = quad(a, b, c)
        # print(f'a: {a}\tb: {b}\tc: {c}\t')
        # print(f't_m: {t_m}\tt_p: {t_p}')

        if t_p > Tube.HIT_THRESHOLD or t_m > Tube.HIT_THRESHOLD:
            # enforce t > 0 (np.nan > 0 == Flase)
            # only time a split (t_m, t_p) -> (+, -) or (-, +)
            # is if ray originates from inside tube and that's not this simulation
            # or reflects of tube
            # so enforce shortest positive distance
            t = np.min([
                t_p if t_p > 0 else np.inf, 
                t_m if t_m > 0 else np.inf, 
            ])
            hit_point = ray(t)

            if np.abs(hit_point[1]) <= self.L / 2:
                return (True, hit_point, t)

        return (False, np.array([np.nan, np.nan, np.nan]), 0)
    
    def reflect(self, ray: Ray, point: np.array) -> Ray:
        # print('TUBE REFLECT')
        point_xz = np.array([point[0], point[2]])
        n_hat = point_xz / linalg.norm(point_xz)
        new_r_hat_xz = reflect_2d(n_hat, ray.hat)
        new = np.array([new_r_hat_xz[0], ray.hat[1], new_r_hat_xz[1]])
        r_hat = new / linalg.norm(new)
        return Ray(ray.lam, point, r_hat)

#####################
# emitting geometry #
#####################

class Sun:
    T    = 5777 # [K]
    flux = 1366 # [W/m^2]

    def __init__(self, solar_zenith: float, trough: Trough, L_air: float = L_AIR):
        """rectangle with the same size as trough to emit from
        
        args:  
            solar_zenith - float: zenith angle of sun [rad]
            trough       - Trough: the current trough
            L_air        - float: modeled length of air at zenith, default L_AIR = 50 [km]
        """

        self.solar_zenith = solar_zenith
        self.solar_altitude = np.pi / 2 - solar_zenith

        # sph2rec has 0 azimuth in x+ direction, so south is pi / 2
        self.sun_hat = sph2rec(1, np.pi / 2, self.solar_zenith)

        self.L_air_actual = rayleigh_airmass(self.solar_altitude, L_air)

        # account for depth of parabola
        extra_L = trough.H_t / np.tan(self.solar_altitude)

        box = np.array([
            [ trough.W/2,  trough.L/2 + extra_L, trough.z_t],
            [-trough.W/2,  trough.L/2 + extra_L, trough.z_t],
            [-trough.W/2, -trough.L/2,           trough.z_t],
            [ trough.W/2, -trough.L/2,           trough.z_t],
        ])
        self.emitting_area = trough.W * (trough.L + extra_L)

        self.sun_corners = box + self.L_air_actual * np.c_[[self.sun_hat for _ in range(4)]]
        # print(self.sun_corners)
        self.blackbody = RandomBlackBody(Sun.T, 1, '.\\table.csv')

    def __repr__(self):
        return f'L`_air:\t\t{self.L_air_actual}\nsun_hat:\t{self.sun_hat}\nsun corners:\n{self.sun_corners}'
    
    def generate_ray(self):
        x = random_range(self.sun_corners[1,0], self.sun_corners[0,0])
        y = random_range(self.sun_corners[2,1], self.sun_corners[0,1])

        return Ray(
            # random wavelength from blackbody table
            self.blackbody.get_lambda(),
            # random x position, random y position + the distance through the air in y, z position of sun
            np.array([x, y, self.sun_corners[1,2]]),
            # ray travels from the sun's direction
            -self.sun_hat
        )