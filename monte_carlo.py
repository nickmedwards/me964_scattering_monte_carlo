from geometry import *
from mc_utils import *

def hit_trace(_ray: Ray, sky: Sky, trough: Trough, tube: Tube):
    trace = [
        _ray.lam, # lambda [mum]
        0,        # absorbed by tube         -> 1
        0,        # absorbed by trough       -> 2
        0,        # absorbed by ground       -> 3
        0,        # absorbed by sky          -> 4
        0,        # tube reflection events   -> 5
        0,        # trough reflection events -> 6
        0,        # scattering events        -> 7
        0,        # absorption events        -> 8
    ]
    # print(np.array(trace[1:5]).sum() > 0)
    ray = _ray
    while np.array(trace[1:5]).sum() == 0:
        # if the ray is pointing down, get ground length, if pointing up get sky length
        # print(ray)
        l_end = get_l_ground(ray) if ray.hat[2] < 0 else sky.get_l_sky(ray)
        l_s = random_l_s(ray.lam)
        l_a = random_l_a(ray.lam)
        tube_hit = tube.is_hit(ray)
        trough_hit = trough.is_hit(ray)
        # print(f'lam:\t{ray.lam}')
        # print(f'l_s:\t{l_s}')
        # print(f'l_a:\t{l_a}')
        # print(f'l_end:\t{l_end}\t{ray.hat[2]}')
        # print(f'tube:\t{tube_hit}')
        # print(f'trough:\t{trough_hit}')
        # print(f'hit selection: {np.argmin([l_end, l_s, l_a, tube_hit[2] if tube_hit[0] else np.inf, trough_hit[2] if trough_hit[0] else np.inf])}')
        min_idx = np.argmin([
            l_end, 
            l_s, 
            l_a,
            tube_hit[2] if tube_hit[0] else np.inf,
            trough_hit[2] if trough_hit[0] else np.inf
        ])
        # if ray.origin[2] < GROUND:
        #     print(f'min_idx: {min_idx}')

        if min_idx == 0: # ground or sky's edge is closest
            trace[3 if ray.hat[2] < 0 else 4] = 1
        elif min_idx == 1: # scattering length is closest
            trace[7] += 1
            # point_s = ray(l_s)
            # psi_s = random_phi()
            # theta_s = random_theta_rayleigh()
            ray = Ray(ray.lam, ray(l_s), sph2rec(1, random_phi(), random_theta_rayleigh()))
        elif min_idx == 2: # absorption length is closest
            trace[8] += 1
            # point_a = ray(l_a)
            # psi_a = random_phi()
            # theta_a = random_theta_emit()
            ray = Ray(ray.lam, ray(l_a), sph2rec(1, random_phi(), random_theta_emit()))
        elif min_idx == 3: # tube is closest
            if is_absorbed(tube.alpha):
                trace[1] = 1
            else:
                ray = tube.reflect(ray, tube_hit[1])
                trace[5] += 1
        elif min_idx == 4: # trough is closest
            if is_reflected(trough.rho):
                ray = trough.reflect(ray, trough_hit[1])
                trace[6] += 1
            else:
                trace[2] = 1

    return trace

def monte_carlo(args):
    func_start = fmt_now()
    N, zenith, m, rho, alpha = args

    trough = Trough(W_T, L_T, m, rho)
    tube = Tube(L_T, 1, alpha)
    sky = Sky()
    sun = Sun(zenith, trough)

    lambdas     = np.zeros(N)
    abs_tube    = np.zeros(N)
    abs_trough  = np.zeros(N)
    abs_ground  = np.zeros(N)
    abs_sky     = np.zeros(N)
    ref_tube    = np.zeros(N)
    ref_trough  = np.zeros(N)
    scatters    = np.zeros(N)
    absorptions = np.zeros(N)

    # energy_ticker = 0
    # time_ticker = 0

    for i in range(N):
        lambdas[i], \
        abs_tube[i], \
        abs_trough[i], \
        abs_ground[i], \
        abs_sky[i], \
        ref_tube[i], \
        ref_trough[i], \
        scatters[i], \
        absorptions[i] = hit_trace(sun.generate_ray(), sky, trough, tube)

        # energy_ticker += (1000 * c_0 * h) / (n_air * lambdas[i])

    df = pd.DataFrame(np.c_[
        lambdas, abs_tube, abs_trough, abs_ground, abs_sky, ref_tube, ref_trough, scatters, absorptions
    ], columns=[
        'lambdas', 'abs_tube', 'abs_trough', 'abs_ground', 'abs_sky', 'ref_tube', 'ref_trough', 'scatters', 'absorptions'
    ])

    # calculate energy of each ray
    # e = c_0 * h / n * lambda, 1000 for converting mum to m
    df['energy'] = (1000 * c_0 * h) / (n_air * df['lambdas']) # [J]
    E_TOTAL  = df['energy'].sum()

    # solar flux [J / (m^2 * s)]
    TIME = E_TOTAL / (Sun.flux * sun.emitting_area) # [s]
    N_ABS_TUBE       = df["abs_tube"].sum() / N
    N_ABS_TROUGH     = df["abs_trough"].sum() / N
    N_ABS_GROUND     = df["abs_ground"].sum() / N
    N_ABS_SKY        = df["abs_sky"].sum() / N
    N_REF_TUBE       = df["ref_tube"].sum()
    N_REF_TROUGH     = df["ref_trough"].sum()
    N_SCATTERS       = df["scatters"].sum()
    mean_SCATTERS    = df["scatters"].mean()
    std_SCATTERS     = df["scatters"].std()
    N_ABSORPTIONS    = df["absorptions"].sum()
    mean_ABSORPTIONS = df["absorptions"].mean()
    std_ABSORPTIONS  = df["absorptions"].std()
    E_TUBE           = df['energy'][df["abs_tube"] == 1].sum()      # [J]
    flux_TUBE        = E_TUBE / (TIME * tube.surface_area)          # [W/m^2]
    mean_lambda_TUBE = df['lambdas'][df["abs_tube"] == 1].mean()
    std_lambda_TUBE  = df['lambdas'][df["abs_tube"] == 1].std()
    E_TROUGH         = df['energy'][df["abs_trough"] == 1].sum()    # [J]
    flux_TROUGH      = E_TROUGH / (TIME * trough.surface_area) # [W/m^2]
    E_GROUND         = df['energy'][df["abs_ground"] == 1].sum()    # [J]
    E_SKY            = df['energy'][df["abs_sky"] == 1].sum()       # [J]
    mean_lambda_SKY  = df['lambdas'][df["abs_sky"] == 1].mean()
    std_lambda_SKY   = df['lambdas'][df["abs_sky"] == 1].std()

    return func_start, np.array([
        N_ABS_TUBE,
        N_ABS_TROUGH,
        N_ABS_GROUND,
        N_ABS_SKY,
        N_REF_TUBE,
        N_REF_TROUGH,
        N_SCATTERS,
        mean_SCATTERS,
        std_SCATTERS,
        N_ABSORPTIONS,
        mean_ABSORPTIONS,
        std_ABSORPTIONS,
        E_TOTAL,
        E_TUBE,
        flux_TUBE,
        mean_lambda_TUBE,
        std_lambda_TUBE,
        E_TROUGH,
        flux_TROUGH,
        E_GROUND,
        E_SKY,
        mean_lambda_SKY,
        std_lambda_SKY,
        TIME
    ])