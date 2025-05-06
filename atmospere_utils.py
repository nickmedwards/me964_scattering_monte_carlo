from math_utils import *
import pvlib

##########################
# set up solana location #
##########################

# lat, long = 32.9235401, -112.9807829
# long_st   = -105 # MST standard meridian
# year      = 2024

# solana = pvlib.location.Location(lat, long)
# solana_area = 7111700 # [m^2]

# zero_pad = lambda n: f'{(3 - len(str(n))) * "0"}{n}'

# srt = solana.get_sun_rise_set_transit(
#     pd.DatetimeIndex(
#         [
#             (datetime.strptime(f'2025 12 {zero_pad(i+1)}', '%Y %H %j').astimezone(timezone(timedelta(hours=-7), name='MST'))).isoformat()
#             for i in range(365)
#         ],
#         tz="MST"
#     )
# )

# monthly_temp = [
#     15, 18, 22, 26, 30, 36, 37, 36, 34, 28, 21, 15
# ]

# # apparent zenith for each day of the year
# apparent_zenith = np.array([
#     solana.get_solarposition(
#         [row['transit'].isoformat()], temperature=monthly_temp[row['transit'].month - 1]
#     )['apparent_zenith']
#     for i, row in srt.iterrows()
# ])


# from apparent_zenith.min()
MIN_ZENITH = 9.483353145565061

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

"""
[1] R. Matichuk et al., “A Decade of Aerosol and Gas Precursor Chemical 
    Characterization at Mt. Lemmon, Arizona (1992 to 2002),” Journal of 
    the Meteorological Society of Japan, vol. 84, no. 4, pp. 653–670, 
    2006, doi: 10.2151/jmsj.84.653.
"""
# single scattering coeffiecient
omega = .937 # at .5 [mum]

"""
[2] H. Horvath, “Spectral extinction coefficients of background aerosols
    in Europe, North and South America: A comparison,” Atmospheric 
    Environment. Part A. General Topics, vol. 25, no. 3–4, pp. 725–732, 
    Jan. 1991, doi: 10.1016/0960-1686(91)90071-E.
"""
# wavelengths used in extinction coefficient measurements
ext_lambda = np.array([.4, .5, .6, .7, .8]) # [mum]
# measured extinction coefficients
beta_km_lambda_measured = np.array([.04, .023, .014, .0085, .0055]) # [km^-1], measured
beta_km_lambda_particles = np.array([.0065, .0084, .0076, .005, .0035]) # [km^-1], smog particles
beta_km_lambda_rayleigh = beta_km_lambda_measured - beta_km_lambda_particles # [km^-1], rayleigh extintion coefficient
beta_choices = [beta_km_lambda_measured, beta_km_lambda_particles, beta_km_lambda_rayleigh]
beta_labels = [('measured', 'Measured'), ('particles', 'Smog Particles'), ('rayleigh', 'Rayleigh Scattering')]
USE_BETA = 0
beta_km_lambda = beta_choices[USE_BETA]

beta_lambda = beta_km_lambda / 1000 # [m^-1]

# measured scattering and absoprtion coefficients as a function of lambda
sigma_s_lambda   = beta_lambda * omega
kappa_lambda     = beta_lambda * (1 - omega)

# added 0 and infinity to avoid bounds issues with f_lerp
x_ext_lambda     = xify(ext_lambda)
f_beta_lambda    = functionify(beta_lambda)
f_sigma_s_lambda = functionify(sigma_s_lambda)
f_kappa_lambda   = functionify(kappa_lambda)

get_beta    = lambda _lambda: f_lerp(_lambda, x_ext_lambda, f_beta_lambda)
get_sigma_s = lambda _lambda: f_lerp(_lambda, x_ext_lambda, f_sigma_s_lambda)
get_kappa   = lambda _lambda: f_lerp(_lambda, x_ext_lambda, f_kappa_lambda)

is_absorbed  = lambda alpha: np.random.random() < alpha
is_reflected = lambda rho: np.random.random() < rho
