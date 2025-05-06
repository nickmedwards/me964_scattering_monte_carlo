from monte_carlo import *

#########################
# convergence study     #
# takes roughly 10 mins #
# use N_RAYS = 10000    #
#########################
if __name__ == '__main__':
    test_trough = Trough(W_T, L_T, M_PARABOLA[1], RHO_PARABOLA[1])
    tube = Tube(L_T, 1, ALPHA_TUBE[1])
    sky = Sky()
    sun = Sun(np.pi / 8, test_trough)

    converge_n_rays = [100, 1000, 10000, 100000]

    CONV_N_ABS_TUBE   = np.zeros((4, 25))
    CONV_N_ABS_TROUGH = np.zeros((4, 25))
    CONV_N_ABS_GROUND = np.zeros((4, 25))
    CONV_N_ABS_SKY    = np.zeros((4, 25))
    CONV_E_TUBE       = np.zeros((4, 25))
    CONV_E_TROUGH     = np.zeros((4, 25))
    CONV_E_GROUND     = np.zeros((4, 25))
    CONV_E_SKY        = np.zeros((4, 25))
    CONV_E_TOTAL      = np.zeros((4, 25))

    for i, n in enumerate(converge_n_rays):
        print(f'N_RAYS  : {n} : {fmt_now()}')
        for j in range(25):
            if n == converge_n_rays[-1]: print(f'100,000 : {j + 1} : {fmt_now()}')
            lambdas     = np.zeros(n)
            abs_tube    = np.zeros(n)
            abs_trough  = np.zeros(n)
            abs_ground  = np.zeros(n)
            abs_sky     = np.zeros(n)
            ref_tube    = np.zeros(n)
            ref_trough  = np.zeros(n)
            scatters    = np.zeros(n)
            absorptions = np.zeros(n)

            for _k in range(n):
                lambdas[_k], \
                abs_tube[_k], \
                abs_trough[_k], \
                abs_ground[_k], \
                abs_sky[_k], \
                ref_tube[_k], \
                ref_trough[_k], \
                scatters[_k], \
                absorptions[_k] = hit_trace(sun.generate_ray(), sky, test_trough, tube)

            df = pd.DataFrame(np.c_[
                lambdas, abs_tube, abs_trough, abs_ground, abs_sky, ref_tube, ref_trough, scatters, absorptions
            ], columns=[
                'lambdas', 'abs_tube', 'abs_trough', 'abs_ground', 'abs_sky', 'ref_tube', 'ref_trough', 'scatters', 'absorptions'
            ])

            df['energy'] = (1000 * c_0 * h) / (n_air * df['lambdas']) # [J]
            e_tot  = df['energy'].sum()

            CONV_N_ABS_TUBE[i, j]    = df["abs_tube"].sum() / n
            CONV_N_ABS_TROUGH[i, j]  = df["abs_trough"].sum() / n
            CONV_N_ABS_GROUND[i, j]  = df["abs_ground"].sum() / n
            CONV_N_ABS_SKY[i, j]     = df["abs_sky"].sum() / n

            CONV_E_TUBE[i, j]   = df['energy'][df["abs_tube"] == 1].sum() / e_tot
            CONV_E_TROUGH[i, j] = df['energy'][df["abs_trough"] == 1].sum() / e_tot
            CONV_E_GROUND[i, j] = df['energy'][df["abs_ground"] == 1].sum() / e_tot
            CONV_E_SKY[i, j]    = df['energy'][df["abs_sky"] == 1].sum() / e_tot
            CONV_E_TOTAL[i, j]  = e_tot / n

    mean_N_ABS_TUBE = np.mean(CONV_N_ABS_TUBE, axis=1)
    mean_N_ABS_TROUGH = np.mean(CONV_N_ABS_TROUGH, axis=1)
    mean_N_ABS_GROUND = np.mean(CONV_N_ABS_GROUND, axis=1)
    mean_N_ABS_SKY = np.mean(CONV_N_ABS_SKY, axis=1)
    mean_E_TUBE = np.mean(CONV_E_TUBE, axis=1)
    mean_E_TROUGH = np.mean(CONV_E_TROUGH, axis=1)
    mean_E_GROUND = np.mean(CONV_E_GROUND, axis=1)
    mean_E_SKY = np.mean(CONV_E_SKY, axis=1)
    mean_E_TOTAL = np.mean(CONV_E_TOTAL, axis=1)

    std_N_ABS_TUBE = np.std(CONV_N_ABS_TUBE, axis=1)
    std_N_ABS_TROUGH = np.std(CONV_N_ABS_TROUGH, axis=1)
    std_N_ABS_GROUND = np.std(CONV_N_ABS_GROUND, axis=1)
    std_N_ABS_SKY = np.std(CONV_N_ABS_SKY, axis=1)
    std_E_TUBE = np.std(CONV_E_TUBE, axis=1)
    std_E_TROUGH = np.std(CONV_E_TROUGH, axis=1)
    std_E_GROUND = np.std(CONV_E_GROUND, axis=1)
    std_E_SKY = np.std(CONV_E_SKY, axis=1)
    std_E_TOTAL = np.std(CONV_E_TOTAL, axis=1)

    with open('.\convergence.txt',mode='w') as f:
        f.writelines([f'{s} : {v}\n' for s, v in [
            ('mean_N_ABS_TUBE', mean_N_ABS_TUBE),
            ('mean_N_ABS_TROUGH', mean_N_ABS_TROUGH),
            ('mean_N_ABS_GROUND', mean_N_ABS_GROUND),
            ('mean_N_ABS_SKY', mean_N_ABS_SKY),
            ('mean_E_TUBE', mean_E_TUBE),
            ('mean_E_TROUGH', mean_E_TROUGH),
            ('mean_E_GROUND', mean_E_GROUND),
            ('mean_E_SKY', mean_E_SKY),
            ('mean_E_TOTAL', mean_E_TOTAL),
            ('std_N_ABS_TUBE', std_N_ABS_TUBE),
            ('std_N_ABS_TROUGH', std_N_ABS_TROUGH),
            ('std_N_ABS_GROUND', std_N_ABS_GROUND),
            ('std_N_ABS_SKY', std_N_ABS_SKY),
            ('std_E_TUBE', std_E_TUBE),
            ('std_E_TROUGH', std_E_TROUGH),
            ('std_E_GROUND', std_E_GROUND),
            ('std_E_SKY', std_E_SKY),
            ('std_E_TOTAL', std_E_TOTAL)
        ]])
    # print(f'N_ABS_TUBE: {CONV_N_ABS_TUBE}')
    # print(f'mean_N_ABS_TUBE: {mean_N_ABS_TUBE}')
    # print(f'std_N_ABS_TUBE: {std_N_ABS_TUBE}')
    # print(f'N_ABS_TROUGH: {CONV_N_ABS_TROUGH}')
    # print(f'mean_N_ABS_TROUGH: {mean_N_ABS_TROUGH}')
    # print(f'std_N_ABS_TROUGH: {std_N_ABS_TROUGH}')
    # print(f'N_ABS_GROUND: {CONV_N_ABS_GROUND}')
    # print(f'mean_N_ABS_GROUND: {mean_N_ABS_GROUND}')
    # print(f'std_N_ABS_GROUND: {std_N_ABS_GROUND}')
    # print(f'N_ABS_SKY: {CONV_N_ABS_SKY}')
    # print(f'mean_N_ABS_SKY: {mean_N_ABS_SKY}')
    # print(f'std_N_ABS_SKY: {std_N_ABS_SKY}')
    # print(f'E_TUBE: {CONV_E_TUBE}')
    # print(f'mean_E_TUBE: {mean_E_TUBE}')
    # print(f'std_E_TUBE: {std_E_TUBE}')
    # print(f'E_TROUGH: {CONV_E_TROUGH}')
    # print(f'mean_E_TROUGH: {mean_E_TROUGH}')
    # print(f'std_E_TROUGH: {std_E_TROUGH}')
    # print(f'E_GROUND: {CONV_E_GROUND}')
    # print(f'mean_E_GROUND: {mean_E_GROUND}')
    # print(f'std_E_GROUND: {std_E_GROUND}')
    # print(f'E_SKY: {CONV_E_SKY}')
    # print(f'mean_E_SKY: {mean_E_SKY}')
    # print(f'std_E_SKY: {std_E_SKY}')
    # print(f'E_TOTAL: {CONV_E_TOTAL}')
    # print(f'mean_E_TOTAL: {mean_E_TOTAL}')
    # print(f'std_E_TOTAL: {std_E_TOTAL}')