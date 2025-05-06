import cProfile, pstats, io
from pstats import SortKey
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from monte_carlo import *

mpl.rcParams['font.size'] = 22

zeniths = [np.pi / 8, np.pi / 4, 3 * np.pi / 8]
results = {}

for zenith in zeniths:
    test_trough = Trough(W_T, L_T, M_PARABOLA[1], RHO_PARABOLA[1])
    tube = Tube(L_T, 1, ALPHA_TUBE[1])
    sky = Sky()
    sun = Sun(zenith, test_trough)

    ###############################
    # profile hit_trace algorithm #
    ###############################


    lambdas     = np.zeros(N_RAYS)
    abs_tube    = np.zeros(N_RAYS)
    abs_trough  = np.zeros(N_RAYS)
    abs_ground  = np.zeros(N_RAYS)
    abs_sky     = np.zeros(N_RAYS)
    ref_tube    = np.zeros(N_RAYS)
    ref_trough  = np.zeros(N_RAYS)
    scatters    = np.zeros(N_RAYS)
    absorptions = np.zeros(N_RAYS)
    # with cProfile.Profile() as pr:
    for i in range(N_RAYS):
        lambdas[i], \
        abs_tube[i], \
        abs_trough[i], \
        abs_ground[i], \
        abs_sky[i], \
        ref_tube[i], \
        ref_trough[i], \
        scatters[i], \
        absorptions[i] = hit_trace(sun.generate_ray(), sky, test_trough, tube)

    df = pd.DataFrame(np.c_[
        lambdas, abs_tube, abs_trough, abs_ground, abs_sky, ref_tube, ref_trough, scatters, absorptions
    ], columns=[
        'lambdas', 'abs_tube', 'abs_trough', 'abs_ground', 'abs_sky', 'ref_tube', 'ref_trough', 'scatters', 'absorptions'
    ])

    # # pr.print_stats()
    # s = io.StringIO()
    # sortby = SortKey.CUMULATIVE
    # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # ps.print_stats()
    # print(s.getvalue())

    # # calculate energy of each ray
    # # e = c_0 * h / n * lambda, 1000 for converting mum to m
    df['energy'] = (1000 * c_0 * h) / (n_air * df['lambdas']) # [J]
    results[zenith] = df
    # E_TOTAL  = df['energy'].sum()

    # # solar flux [J / (m^2 * s)]
    # TIME = E_TOTAL / (Sun.flux * sun.emitting_area) # [s]
    # N_ABS_TUBE       = df["abs_tube"].sum() / N_RAYS
    # N_ABS_TROUGH     = df["abs_trough"].sum() / N_RAYS
    # N_ABS_GROUND     = df["abs_ground"].sum() / N_RAYS
    # N_ABS_SKY        = df["abs_sky"].sum() / N_RAYS
    # N_REF_TUBE       = df["ref_tube"].sum()
    # N_REF_TROUGH     = df["ref_trough"].sum()
    # N_SCATTERS       = df["scatters"].sum()
    # mean_SCATTERS    = df["scatters"].mean()
    # std_SCATTERS     = df["scatters"].std()
    # N_ABSORPTIONS    = df["absorptions"].sum()
    # mean_ABSORPTIONS = df["absorptions"].mean()
    # std_ABSORPTIONS  = df["absorptions"].std()
    # E_TUBE           = df['energy'][df["abs_tube"] == 1].sum()      # [J]
    # flux_TUBE        = E_TUBE / (TIME * tube.surface_area)          # [W/m^2]
    # mean_lambda_TUBE = df['lambdas'][df["abs_tube"] == 1].mean()
    # std_lambda_TUBE  = df['lambdas'][df["abs_tube"] == 1].std()
    # E_TROUGH         = df['energy'][df["abs_trough"] == 1].sum()    # [J]
    # flux_TROUGH      = E_TROUGH / (TIME * test_trough.surface_area) # [W/m^2]
    # E_GROUND         = df['energy'][df["abs_ground"] == 1].sum()    # [J]
    # E_SKY            = df['energy'][df["abs_sky"] == 1].sum()       # [J]
    # mean_lambda_SKY  = df['lambdas'][df["abs_sky"] == 1].mean()
    # std_lambda_SKY   = df['lambdas'][df["abs_sky"] == 1].std()

    # print(f'\n\n{df}\n\n')
    # print(f'time:\t\t{TIME}\t[s]')
    # print(f'N_ABS_TUBE:\t{N_ABS_TUBE}')
    # print(f'N_ABS_TROUGH:\t{N_ABS_TROUGH}')
    # print(f'N_ABS_GROUND:\t{N_ABS_GROUND}')
    # print(f'N_ABS_SKY:\t{N_ABS_SKY}')
    # print(f'N_REF_TUBE:\t{N_REF_TUBE}')
    # print(f'N_REF_TROUGH:\t{N_REF_TROUGH}')
    # print(f'N_SCATTERS:\t{N_SCATTERS}')
    # print(f'mean_SCATTERS:\t{mean_SCATTERS}')
    # print(f'N_ABSORPTIONS:\t{N_ABSORPTIONS}')
    # print(f'E_TUBE:\t\t{E_TUBE}\t[%]')
    # print(f'flux_TUBE:\t{flux_TUBE}\t[W/m^2]')
    # print(f'E_TROUGH:\t{E_TROUGH}\t[%]')
    # print(f'flux_TROUGH:\t{flux_TROUGH}\t[W/m^2]')
    # print(f'E_GROUND:\t{E_GROUND}\t[%]')
    # print(f'E_SKY:\t\t{E_SKY}\t[%]')
    # print(f'ENERGY:\t\t{E_TOTAL}\t[J]')


#########
# plots #
#########

fig, ax_lambda_emitted = plt.subplots(1,1)
fig, ax_scatter_dist = plt.subplots(1,1)
fig, ax_abs_energy_tube = plt.subplots(1,1)
fig, ax_abs_energy_ground = plt.subplots(1,1)
fig, ax_abs_energy_sky = plt.subplots(1,1)
# fig, ax_abs_energy_trough = plt.subplots(1,1)
fig, ax_abs_lambda_tube = plt.subplots(1,1)
fig, ax_abs_lambda_ground = plt.subplots(1,1)
fig, ax_abs_lambda_sky = plt.subplots(1,1)
# fig, ax_abs_lambda_trough = plt.subplots(1,1)

for zenith in zeniths[::-1]:
# run profile to get dfs
    _df = results[zenith]
    _label = f'{np.round(r2d * zenith, 1)} ' + r'[$^{\circ}$]'

    ax_lambda_emitted.hist(_df['lambdas'], histtype="step", bins=200, label=_label) # , alpha=0.6

    ax_scatter_dist.hist(_df['scatters'], histtype="step", bins=int(_df['scatters'].max()), align='left', label=_label) # , alpha=0.6
    ax_scatter_dist.set_yscale('log')
    ax_scatter_dist.set_xbound(-.5, 16.5)

    ax_abs_energy_tube.hist(_df[_df['abs_tube'] == 1]['energy'], histtype="step", bins=100, label=_label) # , alpha=0.6
    ax_abs_energy_ground.hist(_df[_df['abs_ground'] == 1]['energy'], histtype="step", bins=100, label=_label) # , alpha=0.6
    ax_abs_energy_sky.hist(_df[_df['abs_sky'] == 1]['energy'], histtype="step", bins=100, label=_label) # , alpha=0.6
    # ax_abs_energy_trough.hist(_df[_df['abs_trough'] == 1]['energy'], histtype="step", bins=100, label=_label) # , alpha=0.6

    ax_abs_lambda_tube.hist(_df[_df['abs_tube'] == 1]['lambdas'], histtype="step", bins=200, label=_label) # , alpha=0.6
    ax_abs_lambda_ground.hist(_df[_df['abs_ground'] == 1]['lambdas'], histtype="step", bins=200, label=_label) # , alpha=0.6
    ax_abs_lambda_sky.hist(_df[_df['abs_sky'] == 1]['lambdas'], histtype="step", bins=200, label=_label) # , alpha=0.6
    # ax_abs_lambda_trough.hist(_df[_df['abs_trough'] == 1]['lambdas'], histtype="step", bins=200, label=_label) # , alpha=0.6



ax_lambda_emitted.set_title('Wavelength Distribution Emitted\nat Selected Zenith Angles')
ax_lambda_emitted.set_xlabel(r'Wavelength [$\mu$m]')
ax_lambda_emitted.set_ylabel('Frequency')
ax_lambda_emitted.legend()

ax_scatter_dist.set_title('Distribution of Scattering Events per Ray\nat Selected Zenith Angles')
ax_scatter_dist.set_xlabel('Scattering Events [-]')
ax_scatter_dist.set_ylabel('Frequency')
ax_scatter_dist.legend()

ax_abs_energy_tube.set_title('Energy Distribution Absorbed by Tube\nat Selected Zenith Angles')
ax_abs_energy_tube.set_xlabel('Energy [J]')
ax_abs_energy_tube.set_ylabel('Frequency')
ax_abs_energy_tube.legend()

ax_abs_energy_ground.set_title('Energy Distribution Absorbed by Ground\nat Selected Zenith Angles')
ax_abs_energy_ground.set_xlabel('Energy [J]')
ax_abs_energy_ground.set_ylabel('Frequency')
ax_abs_energy_ground.legend()

ax_abs_energy_sky.set_title('Energy Distribution Exiting Sky\nat Selected Zenith Angles')
ax_abs_energy_sky.set_xlabel('Energy [J]')
ax_abs_energy_sky.set_ylabel('Frequency')
ax_abs_energy_sky.legend()

# ax_abs_energy_trough.set_title('Energy Distribution Absorbed by Trough\nat Selected Zenith Angles')
# ax_abs_energy_trough.set_xlabel('Energy [J]')
# ax_abs_energy_trough.set_ylabel('Frequency')
# ax_abs_energy_trough.legend()

ax_abs_lambda_tube.set_title('Wavelength Distribution Absorbed by Tube\nat Selected Zenith Angles')
ax_abs_lambda_tube.set_xlabel(r'Wavelength [$\mu$m]')
ax_abs_lambda_tube.set_ylabel('Frequency')
ax_abs_lambda_tube.legend()
ax_abs_lambda_tube.set_xbound(0, 4)

ax_abs_lambda_ground.set_title('Wavelength Distribution Absorbed by Ground\nat Selected Zenith Angles')
ax_abs_lambda_ground.set_xlabel(r'Wavelength [$\mu$m]')
ax_abs_lambda_ground.set_ylabel('Frequency')
ax_abs_lambda_ground.legend()
ax_abs_lambda_ground.set_xbound(0, 4)

ax_abs_lambda_sky.set_title('Wavelength Distribution Exiting Sky\nat Selected Zenith Angles')
ax_abs_lambda_sky.set_xlabel(r'Wavelength [$\mu$m]')
ax_abs_lambda_sky.set_ylabel('Frequency')
ax_abs_lambda_sky.legend()
ax_abs_lambda_sky.set_xbound(0, 4)

# ax_abs_lambda_trough.set_title('Wavelength Distribution Absorbed by Trough\nat Selected Zenith Angles')
# ax_abs_lambda_trough.set_xlabel(r'Wavelength [$\mu$m]')
# ax_abs_lambda_trough.set_ylabel('Frequency')
# ax_abs_lambda_trough.legend()


# ax.set_xbound(0, 4)

# fig, ax = plt.subplots(1,1)


# ax = plt.figure().add_subplot(projection='3d')
# ax.scatter(*test_trough.for_plot())
# ax.scatter(*tube.for_plot())

##############
# test plots #
##############

# # test plot for black body dist
# fig, ax = plt.subplots(1,1)
# T_sun = 5777
# random_emitter = RandomBlackBody(T_sun, 1, '.\\table.csv')
# ax.hist([random_emitter.get_lambda() for _ in range(10000)], bins=200)

# test for extinction fit
# fig, ax = plt.subplots(1,1)
# test_ext_lambda = np.linspace(ext_lambda[0], ext_lambda[-1], num=30)
# ax.plot(test_ext_lambda, [get_sigma_s(l) for l in test_ext_lambda])
# ax.plot(test_ext_lambda, [get_kappa(l) for l in test_ext_lambda])
# ax.scatter(ext_lambda, beta_lambda, c='g')

# test for rayleigh inversion
# fig, ax = plt.subplots(1,1)
# ax.scatter(R_theta_rayleigh, theta_rayleigh, c='g', marker='.')

# trough and tube
# ax = plt.figure().add_subplot(projection='3d')
# ax.scatter(*test_trough.for_plot())
# ax.scatter(*tube.for_plot())

# apparent zenith for noons
# fig, ax = plt.subplots(1,1)
# ax.plot(apparent_zenith)

# sky
# fig, ax = plt.subplots(1,1)
# ax.plot(sky.z, sky.air)

plt.show()