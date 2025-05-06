import concurrent.futures
from monte_carlo import *

if __name__ == '__main__':
    zeniths = d2r * np.linspace(MIN_ZENITH, 89, num=91)
    TEST_MATRIX = np.zeros((len(zeniths) * len(M_PARABOLA) * len(RHO_PARABOLA) * len(ALPHA_TUBE), 5))
    TEST_MATRIX_ARGS = []

    ticker = 0
    for _i in range(len(zeniths)):
        for _j in range(len(M_PARABOLA)):
            for _k in range(len(RHO_PARABOLA)):
                for _l in range(len(ALPHA_TUBE)):
                    TEST_MATRIX[ticker, :] = np.array([
                        N_RAYS, # int(RAYS[_i]),
                        zeniths[_i],
                        M_PARABOLA[_j],
                        RHO_PARABOLA[_k],
                        ALPHA_TUBE[_l]
                    ])
                    TEST_MATRIX_ARGS.append((
                        N_RAYS, # int(RAYS[_i]),
                        zeniths[_i],
                        M_PARABOLA[_j],
                        RHO_PARABOLA[_k],
                        ALPHA_TUBE[_l]
                    ))    
                    ticker += 1

    n_TEST = len(TEST_MATRIX_ARGS)
    arr_N_ABS_TUBE       = np.zeros(n_TEST)
    arr_N_ABS_TROUGH     = np.zeros(n_TEST)
    arr_N_ABS_GROUND     = np.zeros(n_TEST)
    arr_N_ABS_SKY        = np.zeros(n_TEST)
    arr_N_REF_TUBE       = np.zeros(n_TEST)
    arr_N_REF_TROUGH     = np.zeros(n_TEST)

    arr_N_SCATTERS       = np.zeros(n_TEST)
    arr_mean_SCATTERS    = np.zeros(n_TEST)
    arr_std_SCATTERS     = np.zeros(n_TEST)
    
    arr_N_ABSORPTIONS    = np.zeros(n_TEST)
    arr_mean_ABSORPTIONS = np.zeros(n_TEST)
    arr_std_ABSORPTIONS  = np.zeros(n_TEST)
    
    arr_E_TOTAL          = np.zeros(n_TEST)
    
    arr_E_TUBE           = np.zeros(n_TEST)
    arr_flux_TUBE        = np.zeros(n_TEST)
    arr_mean_lambda_TUBE = np.zeros(n_TEST)
    arr_std_lambda_TUBE  = np.zeros(n_TEST)
    
    arr_E_TROUGH         = np.zeros(n_TEST)
    arr_flux_TROUGH      = np.zeros(n_TEST)
    arr_E_GROUND         = np.zeros(n_TEST)
    
    arr_E_SKY            = np.zeros(n_TEST)
    arr_mean_lambda_SKY  = np.zeros(n_TEST)
    arr_std_lambda_SKY   = np.zeros(n_TEST)
    
    arr_TIME             = np.zeros(n_TEST)

    with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
        start = fmt_now()
        for i, vals in enumerate(executor.map(monte_carlo, TEST_MATRIX_ARGS)):
            func_start, arr = vals
            print(f'{i + 1} / {n_TEST} | {np.round(100 * (i + 1) / n_TEST, 4)}\t: {start} -- {func_start} -- {fmt_now()}')
            arr_N_ABS_TUBE[i], \
            arr_N_ABS_TROUGH[i], \
            arr_N_ABS_GROUND[i], \
            arr_N_ABS_SKY[i], \
            arr_N_REF_TUBE[i], \
            arr_N_REF_TROUGH[i], \
            arr_N_SCATTERS[i], \
            arr_mean_SCATTERS[i], \
            arr_std_SCATTERS[i], \
            arr_N_ABSORPTIONS[i], \
            arr_mean_ABSORPTIONS[i], \
            arr_std_ABSORPTIONS[i], \
            arr_E_TOTAL[i], \
            arr_E_TUBE[i], \
            arr_flux_TUBE[i], \
            arr_mean_lambda_TUBE[i], \
            arr_std_lambda_TUBE[i], \
            arr_E_TROUGH[i], \
            arr_flux_TROUGH[i], \
            arr_E_GROUND[i], \
            arr_E_SKY[i], \
            arr_mean_lambda_SKY[i], \
            arr_std_lambda_SKY[i], \
            arr_TIME[i] = arr
        
    df = pd.DataFrame(np.c_[
        *TEST_MATRIX.T,
         arr_N_ABS_TUBE,
         arr_N_ABS_TROUGH,
         arr_N_ABS_GROUND,
         arr_N_ABS_SKY,
         arr_N_REF_TUBE,
         arr_N_REF_TROUGH,
         arr_N_SCATTERS,
         arr_mean_SCATTERS,
         arr_std_SCATTERS,
         arr_N_ABSORPTIONS,
         arr_mean_ABSORPTIONS,
         arr_std_ABSORPTIONS,
         arr_E_TOTAL,
         arr_E_TUBE,
         arr_flux_TUBE,
         arr_mean_lambda_TUBE,
         arr_std_lambda_TUBE,
         arr_E_TROUGH,
         arr_flux_TROUGH,
         arr_E_GROUND,
         arr_E_SKY,
         arr_mean_lambda_SKY,
         arr_std_lambda_SKY,
         arr_TIME
    ], columns=[
        'N_RAYS',
        'zeniths',
        'M_PARABOLA',
        'RHO_PARABOLA',
        'ALPHA_TUBE',
        'arr_N_ABS_TUBE',
        'arr_N_ABS_TROUGH', 
        'arr_N_ABS_GROUND', 
        'arr_N_ABS_SKY', 
        'arr_N_REF_TUBE', 
        'arr_N_REF_TROUGH', 
        'arr_N_SCATTERS', 
        'arr_mean_SCATTERS', 
        'arr_std_SCATTERS',
        'arr_N_ABSORPTIONS', 
        'arr_mean_ABSORPTIONS', 
        'arr_std_ABSORPTIONS', 
        'arr_E_TOTAL', 
        'arr_E_TUBE', 
        'arr_flux_TUBE', 
        'arr_mean_lambda_TUBE', 
        'arr_std_lambda_TUBE', 
        'arr_E_TROUGH', 
        'arr_flux_TROUGH', 
        'arr_E_GROUND', 
        'arr_E_SKY', 
        'arr_mean_lambda_SKY', 
        'arr_std_lambda_SKY', 
        'arr_TIME'
    ])

    df.to_csv(f'.\\parametric_{beta_labels[USE_BETA][0]}.csv')