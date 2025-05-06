import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from geometry import *
from atmospere_utils import beta_labels

mpl.rcParams['font.size'] = 22

# for name in ['measured', 'particles']:
name = 'measured'
plt.style.use('seaborn-v0_8-colorblind')
# rays = '30000'
rays = 'fixed'
# for rays in ['fixed']: # '10000', '30000', 
for (fname, _label) in beta_labels[:1]:
    df = pd.read_csv(f'.\\parametric_{fname}.csv')

    # fig, ax = plt.subplots(1,1)
    # ax.errorbar(df['zeniths'], df['arr_mean_SCATTERS'], df['arr_std_SCATTERS'])
    # ax.set_title(f'{name} - {rays[:2]}k: scatter events')
    # fig, ax = plt.subplots(1,1)
    # ax.scatter(df['zeniths'], df['arr_mean_SCATTERS'])
    # ax.set_title(f'{name} - {rays[:2]}k: scatter events')

    # fig, ax = plt.subplots(1,1)
    # ax.scatter(df[df['ALPHA_TUBE'] == 1.0]['zeniths'], df[df['ALPHA_TUBE'] == 1.0]['arr_TIME'] * df[df['ALPHA_TUBE'] == 1.0]['arr_flux_TUBE'], label='1.0', alpha=0.9)
    # ax.scatter(df[df['ALPHA_TUBE'] == .98]['zeniths'], df[df['ALPHA_TUBE'] == .98]['arr_TIME'] * df[df['ALPHA_TUBE'] == .98]['arr_flux_TUBE'], label='.98', alpha=0.8)
    # ax.scatter(df[df['ALPHA_TUBE'] == .96]['zeniths'], df[df['ALPHA_TUBE'] == .96]['arr_TIME'] * df[df['ALPHA_TUBE'] == .96]['arr_flux_TUBE'], label='.96', alpha=0.7)
    # ax.scatter(df[df['ALPHA_TUBE'] == .94]['zeniths'], df[df['ALPHA_TUBE'] == .94]['arr_TIME'] * df[df['ALPHA_TUBE'] == .94]['arr_flux_TUBE'], label='.94', alpha=0.6)
    # ax.legend()
    # ax.set_title(f'{name} - {rays[:2]}k: alpha tube - energy')

    # fig, ax = plt.subplots(1,1)
    # ax.scatter(df[df['RHO_PARABOLA'] == 1.0]['zeniths'], df[df['RHO_PARABOLA'] == 1.0]['arr_TIME'] * df[df['RHO_PARABOLA'] == 1.0]['arr_flux_TUBE'], label='1.0', alpha=0.9)
    # ax.scatter(df[df['RHO_PARABOLA'] == .98]['zeniths'], df[df['RHO_PARABOLA'] == .98]['arr_TIME'] * df[df['RHO_PARABOLA'] == .98]['arr_flux_TUBE'], label='.98', alpha=0.8)
    # ax.scatter(df[df['RHO_PARABOLA'] == .96]['zeniths'], df[df['RHO_PARABOLA'] == .96]['arr_TIME'] * df[df['RHO_PARABOLA'] == .96]['arr_flux_TUBE'], label='.96', alpha=0.7)
    # ax.scatter(df[df['RHO_PARABOLA'] == .94]['zeniths'], df[df['RHO_PARABOLA'] == .94]['arr_TIME'] * df[df['RHO_PARABOLA'] == .94]['arr_flux_TUBE'], label='.94', alpha=0.6)
    # ax.legend()
    # ax.set_title(f'{name} - {rays[:2]}k: rho trough - energy')

    fig, ax = plt.subplots(1,1)
    ax.scatter(df[df['M_PARABOLA'] == M_PARABOLA[0]]['zeniths'], df[df['M_PARABOLA'] == M_PARABOLA[0]]['arr_E_TUBE'], alpha=1.0, label=f'{M_PARABOLA[0]}', marker='.')
    ax.scatter(df[df['M_PARABOLA'] == M_PARABOLA[1]]['zeniths'], df[df['M_PARABOLA'] == M_PARABOLA[1]]['arr_E_TUBE'], alpha=0.8, label=f'{M_PARABOLA[1]}', marker='.')
    ax.scatter(df[df['M_PARABOLA'] == M_PARABOLA[2]]['zeniths'], df[df['M_PARABOLA'] == M_PARABOLA[2]]['arr_E_TUBE'], alpha=0.7, label=f'{M_PARABOLA[2]}', marker='.')
    ax.scatter(df[df['M_PARABOLA'] == M_PARABOLA[3]]['zeniths'], df[df['M_PARABOLA'] == M_PARABOLA[3]]['arr_E_TUBE'], alpha=0.7, label=f'{M_PARABOLA[3]}', marker='.')
    ax.set_title(f'Tube Flux vs. Zenith Angle\nColorcoded by Trough Geometry')
    ax.legend()
    ax.grid(visible=True, alpha=.7)
    
    # fig, ax = plt.subplots(1,1)
    # ax.scatter(df[df['ALPHA_TUBE'] == ALPHA_TUBE[0]]['zeniths'], df[df['ALPHA_TUBE'] == ALPHA_TUBE[0]]['arr_E_TUBE'], alpha=1.0, label=f'{ALPHA_TUBE[0]}', marker='.')
    # ax.scatter(df[df['ALPHA_TUBE'] == ALPHA_TUBE[1]]['zeniths'], df[df['ALPHA_TUBE'] == ALPHA_TUBE[1]]['arr_E_TUBE'], alpha=0.8, label=f'{ALPHA_TUBE[1]}', marker='.')
    # ax.scatter(df[df['ALPHA_TUBE'] == ALPHA_TUBE[2]]['zeniths'], df[df['ALPHA_TUBE'] == ALPHA_TUBE[2]]['arr_E_TUBE'], alpha=0.7, label=f'{ALPHA_TUBE[2]}', marker='.')
    # ax.scatter(df[df['ALPHA_TUBE'] == ALPHA_TUBE[3]]['zeniths'], df[df['ALPHA_TUBE'] == ALPHA_TUBE[3]]['arr_E_TUBE'], alpha=0.7, label=f'{ALPHA_TUBE[3]}', marker='.')
    # ax.set_title(f'Tube Flux vs. Zenith Angle\nColorcoded by Tube Absorptance')
    # ax.legend()
    # ax.grid(visible=True, alpha=.7)

    # fig, ax = plt.subplots(1,1)
    # ax.scatter(df[df['M_PARABOLA'] == .06]['zeniths'], df[df['M_PARABOLA'] == .06]['arr_flux_TROUGH'], label='.06')
    # ax.scatter(df[df['M_PARABOLA'] == .033]['zeniths'], df[df['M_PARABOLA'] == .033]['arr_flux_TROUGH'], label='.033')
    # ax.scatter(df[df['M_PARABOLA'] == .025]['zeniths'], df[df['M_PARABOLA'] == .025]['arr_flux_TROUGH'], label='.025')
    # ax.scatter(df[df['M_PARABOLA'] == .0175]['zeniths'], df[df['M_PARABOLA'] == .0175]['arr_flux_TROUGH'], label='.0175')
    # ax.set_title(f'{name} - {rays[:2]}k: m parabola - trough flux')
    # ax.legend()


    # fig, ax = plt.subplots(1,1)
    # ax.scatter(df[df['M_PARABOLA'] == .06]['zeniths'], df[df['M_PARABOLA'] == .06]['arr_TIME'], label='.06')
    # ax.scatter(df[df['M_PARABOLA'] == .033]['zeniths'], df[df['M_PARABOLA'] == .033]['arr_TIME'], label='.033')
    # ax.scatter(df[df['M_PARABOLA'] == .025]['zeniths'], df[df['M_PARABOLA'] == .025]['arr_TIME'], label='.025')
    # ax.scatter(df[df['M_PARABOLA'] == .0175]['zeniths'], df[df['M_PARABOLA'] == .0175]['arr_TIME'], label='.0175')
    # ax.legend()
    # # ax.scatter(df['zeniths'], df['arr_TIME'])
    # ax.set_title(f'{name} - {rays[:2]}k: time')

    # fig, ax = plt.subplots(1,1)
    # ax.scatter(df['zeniths'], df['arr_mean_lambda_TUBE'])
    # ax.scatter(df['zeniths'], df['arr_mean_lambda_SKY'])
    # ax.set_ybound(0, 4)
    # ax.set_title(f'{name} - {rays[:2]}k: mean lambdas')

    # fig, ax = plt.subplots(1,1)
    # ax.errorbar(df['zeniths'], df['arr_mean_lambda_TUBE'], df['arr_std_lambda_TUBE'])
    # ax.set_title(f'{name} - {rays[:2]}k: mean and std lambdas tube')

    # fig, ax = plt.subplots(1,1)
    # ax.errorbar(df['zeniths'], df['arr_mean_lambda_SKY'], df['arr_std_lambda_SKY'])
    # ax.set_title(f'{name} - {rays[:2]}k: mean and std lambdas sky')

    # fig, ax = plt.subplots(1,1)
    # ax.scatter(df['zeniths'], df['arr_mean_lambda_TUBE'])
    # ax.set_title(f'{name} - {rays[:2]}k: mean and std lambdas tube')

    # fig, ax = plt.subplots(1,1)
    # ax.scatter(df['zeniths'], df['arr_mean_lambda_SKY'])
    # ax.set_title(f'{name} - {rays[:2]}k: mean and std lambdas sky')

    # fig, ax = plt.subplots(1,1)
    # ax.scatter(df['zeniths'], df['arr_N_ABS_TUBE'], label='Tube [%]')
    # ax.scatter(df['zeniths'], df['arr_N_ABS_SKY'], label='Ground [%]')
    # ax.scatter(df['zeniths'], df['arr_N_ABS_GROUND'], label='Sky [%]')
    # ax.scatter(df['zeniths'], df['arr_N_ABS_TROUGH'], label='Trough [%]')
    # ax.legend()
    # ax.set_title(f'{name} - {rays[:2]}k: N absorptions')

    # fig, ax = plt.subplots(1,1)
    # ax.scatter(df['zeniths'], df['arr_E_TUBE'] / df['arr_E_TOTAL'], label='Tube [%]')
    # ax.scatter(df['zeniths'], df['arr_E_SKY'] / df['arr_E_TOTAL'], label='Ground [%]')
    # ax.scatter(df['zeniths'], df['arr_E_GROUND'] / df['arr_E_TOTAL'], label='Sky [%]')
    # ax.scatter(df['zeniths'], df['arr_E_TROUGH'] / df['arr_E_TOTAL'], label='Trough [%]')
    # ax.legend()
    # ax.set_title(f'{name} - {rays[:2]}k: energy')

    # fig, ax = plt.subplots(1,1)
    # tmp_zeniths = df[(df['M_PARABOLA'] == M_PARABOLA[0]) & (df['RHO_PARABOLA'] == RHO_PARABOLA[0])]['zeniths']
    # tmp_time = df[(df['M_PARABOLA'] == M_PARABOLA[0]) & (df['RHO_PARABOLA'] == RHO_PARABOLA[0])]['arr_TIME']
    # tmp_E = df[(df['M_PARABOLA'] == M_PARABOLA[0]) & (df['RHO_PARABOLA'] == RHO_PARABOLA[0])]['arr_E_TOTAL']
    # sun_area = np.array([Sun(z, Trough(W_T, L_T, M_PARABOLA[0], RHO_PARABOLA[0])).emitting_area for z in tmp_zeniths])

    # ax.scatter(tmp_zeniths, tmp_E / (tmp_time * sun_area))
    # ax.set_title(f'{name} - {rays[:2]}k: sun area')


# fig, ax_zenith_v_scatters = plt.subplots(1,1)
fig, ax_airmass_v_scatters = plt.subplots(1,1)
# fig, ax_zenith_v_flux_TUBE = plt.subplots(1,1)
fig, ax_airmass_v_flux_TUBE = plt.subplots(1,1)
# fig, ax_zenith_v_mean_lambda_TUBE = plt.subplots(1,1)
fig, ax_airmass_v_mean_lambda_TUBE = plt.subplots(1,1)
# fig, ax_scatters_v_mean_lambda_TUBE = plt.subplots(1,1)
# fig, ax_zenith_v_mean_lambda_SKY = plt.subplots(1,1)
fig, ax_airmass_v_mean_lambda_SKY = plt.subplots(1,1)
# fig, ax_scatters_v_mean_lambda_SKY = plt.subplots(1,1)
# fig, ax_zenith_v_std_lambda_TUBE = plt.subplots(1,1)
fig, ax_airmass_v_std_lambda_TUBE = plt.subplots(1,1)
# fig, ax_scatters_v_std_lambda_TUBE = plt.subplots(1,1)
# fig, ax_zenith_v_std_lambda_SKY = plt.subplots(1,1)
fig, ax_airmass_v_std_lambda_SKY = plt.subplots(1,1)
# fig, ax_scatters_v_std_lambda_SKY = plt.subplots(1,1)

for i, pair in enumerate(beta_labels):
    fname, _label = pair
    df = pd.read_csv(f'.\\parametric_{fname}.csv')
    df['air'] = np.array([rayleigh_airmass(np.pi / 2 - z) for z in df['zeniths']])

    # ax_zenith_v_scatters.scatter(r2d * df['zeniths'], df['arr_mean_SCATTERS'], label=_label, marker='.')
    ax_airmass_v_scatters.scatter(df['air'], df['arr_mean_SCATTERS'], label=_label, marker='.')
    
    # ax_zenith_v_flux_TUBE.scatter(r2d * df['zeniths'], df['arr_flux_TUBE'], label=_label, marker='.')
    ax_airmass_v_flux_TUBE.scatter(df['air'], df['arr_flux_TUBE'], label=_label, marker='.')
    
    # ax_zenith_v_mean_lambda_TUBE.scatter(r2d * df['zeniths'], df['arr_mean_lambda_TUBE'], label=_label, marker='.')
    ax_airmass_v_mean_lambda_TUBE.scatter(df['air'], df['arr_mean_lambda_TUBE'], label=_label, marker='.')
    # ax_scatters_v_mean_lambda_TUBE.scatter(df['arr_mean_SCATTERS'], df['arr_mean_lambda_TUBE'], label=_label, marker='.', alpha=1 if i < 1 else .7)

    # ax_zenith_v_mean_lambda_SKY.scatter(r2d * df['zeniths'], df['arr_mean_lambda_SKY'], label=_label, marker='.')
    ax_airmass_v_mean_lambda_SKY.scatter(df['air'], df['arr_mean_lambda_SKY'], label=_label, marker='.')
    # ax_scatters_v_mean_lambda_SKY.scatter(df['arr_mean_SCATTERS'], df['arr_mean_lambda_SKY'], label=_label, marker='.', alpha=1 if i < 1 else .7)
    
    # ax_zenith_v_std_lambda_TUBE.scatter(r2d * df['zeniths'], df['arr_std_lambda_TUBE'], label=_label, marker='.', alpha=1 if i < 1 else .7)
    ax_airmass_v_std_lambda_TUBE.scatter(df['air'], df['arr_std_lambda_TUBE'], label=_label, marker='.', alpha=1 if i < 1 else .5)
    # ax_scatters_v_std_lambda_TUBE.scatter(df['arr_mean_SCATTERS'], df['arr_std_lambda_TUBE'], label=_label, marker='.', alpha=1 if i < 1 else .7)
    
    # ax_zenith_v_std_lambda_SKY.scatter(r2d * df['zeniths'], df['arr_std_lambda_SKY'], label=_label, marker='.', alpha=1 if i < 1 else .7)
    ax_airmass_v_std_lambda_SKY.scatter(df['air'], df['arr_std_lambda_SKY'], label=_label, marker='.', alpha=1 if i < 1 else .5)
    # ax_scatters_v_std_lambda_SKY.scatter(df['arr_mean_SCATTERS'], df['arr_std_lambda_SKY'], label=_label, marker='.', alpha=1 if i < 1 else .7)

    # ax_std_lambda_TUBE.scatter(r2d * df['zeniths'], df['arr_std_lambda_TUBE'], label=name)
    # ax_std_lambda_SKY.scatter(r2d * df['zeniths'], df['arr_std_lambda_SKY'], label=name)


    # axs[0].scatter(df[df['M_PARABOLA'] == .06]['zeniths'], df[df['M_PARABOLA'] == .06]['arr_TIME'], label=f'{name} - .06')
    # axs[0].scatter(df[df['M_PARABOLA'] == .033]['zeniths'], df[df['M_PARABOLA'] == .033]['arr_TIME'], label=f'{name} - .033')
    # axs[0].scatter(df[df['M_PARABOLA'] == .025]['zeniths'], df[df['M_PARABOLA'] == .025]['arr_TIME'], label=f'{name} - .025')
    # axs[0].scatter(df[df['M_PARABOLA'] == .0175]['zeniths'], df[df['M_PARABOLA'] == .0175]['arr_TIME'], label=f'{name} - .0175')
    # axs[1].scatter(df[df['M_PARABOLA'] == .06]['zeniths'], df[df['M_PARABOLA'] == .06]['arr_flux_TUBE'], alpha=1 if i ==0 else .5, marker='.')
    # axs[1].scatter(df[df['M_PARABOLA'] == .033]['zeniths'], df[df['M_PARABOLA'] == .033]['arr_flux_TUBE'], alpha=1 if i ==0 else .4, marker='.')
    # axs[1].scatter(df[df['M_PARABOLA'] == .025]['zeniths'], df[df['M_PARABOLA'] == .025]['arr_flux_TUBE'], alpha=1 if i ==0 else .3, marker='.')
    # axs[1].scatter(df[df['M_PARABOLA'] == .0175]['zeniths'], df[df['M_PARABOLA'] == .0175]['arr_flux_TUBE'], alpha=1 if i ==0 else .2, marker='.')

deg = r'[$^{\circ}$]'
flux = r'[W/m$^{2}$]'
mum = r'[$\mu$m]'

# ax_zenith_v_scatters.set_title(f'Zenith Angle vs. Scattering Events')
# ax_zenith_v_scatters.set_xlabel(f'Zenith Angle {deg}')
# ax_zenith_v_scatters.set_ylabel(f'Scattering Events [-]')
# ax_zenith_v_scatters.legend()
# ax_zenith_v_scatters.grid(visible=True, alpha=.7)

ax_airmass_v_scatters.set_title(f'Air Length vs. Scattering Events')
ax_airmass_v_scatters.set_xlabel(f'Air Length [m]')
ax_airmass_v_scatters.set_ylabel(f'Scattering Events [-]')
ax_airmass_v_scatters.legend()
ax_airmass_v_scatters.grid(visible=True, alpha=.7)
ax_airmass_v_scatters.set_xscale('log')

# ax_zenith_v_flux_TUBE.set_title(f'Zenith Angle vs. Radiation Flux on Tube')
# ax_zenith_v_flux_TUBE.set_xlabel(f'Zenith Angle {deg}')
# ax_zenith_v_flux_TUBE.set_ylabel(f'Radiation Flux on Tube {flux}')
# ax_zenith_v_flux_TUBE.legend()
# ax_zenith_v_flux_TUBE.grid(visible=True, alpha=.7)

ax_airmass_v_flux_TUBE.set_title(f'Air Length vs. Radiation Flux on Tube')
ax_airmass_v_flux_TUBE.set_xlabel(f'Air Length [m]')
ax_airmass_v_flux_TUBE.set_ylabel(f'Radiation Flux on Tube {flux}')
ax_airmass_v_flux_TUBE.legend()
ax_airmass_v_flux_TUBE.grid(visible=True, alpha=.7)
ax_airmass_v_flux_TUBE.set_xscale('log')

# ax_zenith_v_mean_lambda_TUBE.set_title(f'Zenith Angle vs. Mean Wavelength Absorbed by Tube')
# ax_zenith_v_mean_lambda_TUBE.set_xlabel(f'Zenith Angle {deg}')
# ax_zenith_v_mean_lambda_TUBE.set_ylabel(f'Mean Wavelength Absorbed by Tube {mum}')
# ax_zenith_v_mean_lambda_TUBE.legend()
# ax_zenith_v_mean_lambda_TUBE.grid(visible=True, alpha=.7)
# ax_zenith_v_mean_lambda_TUBE.set_ybound(0.75, 2)

ax_airmass_v_mean_lambda_TUBE.set_title(f'Air Length vs. Mean Wavelength Absorbed by Tube')
ax_airmass_v_mean_lambda_TUBE.set_xlabel(f'Air Length [m]')
ax_airmass_v_mean_lambda_TUBE.set_ylabel(f'Mean Wavelength Absorbed by Tube {mum}')
ax_airmass_v_mean_lambda_TUBE.legend()
ax_airmass_v_mean_lambda_TUBE.grid(visible=True, alpha=.7)
ax_airmass_v_mean_lambda_TUBE.set_ybound(0.75, 2)
ax_airmass_v_mean_lambda_TUBE.set_xscale('log')

# ax_scatters_v_mean_lambda_TUBE.set_title(f'Mean Scattering Events vs.\nMean Wavelgth Absorbed by Tube')
# ax_scatters_v_mean_lambda_TUBE.set_xlabel(f'Mean Scattering Events[-]')
# ax_scatters_v_mean_lambda_TUBE.set_ylabel(f'Mean Wavelgth Absorbed by Tube {mum}')
# ax_scatters_v_mean_lambda_TUBE.legend()
# ax_scatters_v_mean_lambda_TUBE.grid(visible=True, alpha=.7)
# ax_scatters_v_mean_lambda_TUBE.set_ybound(0.75, 2)

# ax_zenith_v_mean_lambda_SKY.set_title(f'Zenith Angle vs. Mean Wavelength Exiting Sky')
# ax_zenith_v_mean_lambda_SKY.set_xlabel(f'Zenith Angle {deg}')
# ax_zenith_v_mean_lambda_SKY.set_ylabel(f'Mean Wavelength Exiting Sky {mum}')
# ax_zenith_v_mean_lambda_SKY.legend()
# ax_zenith_v_mean_lambda_SKY.grid(visible=True, alpha=.7)

ax_airmass_v_mean_lambda_SKY.set_title(f'Air Length vs. Mean Wavelength Exiting Sky')
ax_airmass_v_mean_lambda_SKY.set_xlabel(f'Air Length [m]')
ax_airmass_v_mean_lambda_SKY.set_ylabel(f'Mean Wavelength Exiting Sky {mum}')
ax_airmass_v_mean_lambda_SKY.legend()
ax_airmass_v_mean_lambda_SKY.grid(visible=True, alpha=.7)
ax_airmass_v_mean_lambda_SKY.set_xscale('log')

# ax_scatters_v_mean_lambda_SKY.set_title(f'Mean Scattering Events vs. Mean Wavelgth Absorbed Exiting Sky')
# ax_scatters_v_mean_lambda_SKY.set_xlabel(f'Mean Scattering Events [-]')
# ax_scatters_v_mean_lambda_SKY.set_ylabel(f'Mean Wavelgth Absorbed Exiting Sky {mum}')
# ax_scatters_v_mean_lambda_SKY.legend()
# ax_scatters_v_mean_lambda_SKY.grid(visible=True, alpha=.7)

# ax_zenith_v_std_lambda_TUBE.set_title(f'Zenith Angle vs. Stadard Deviation\nof Wavelength Absorbed by Tube')
# ax_zenith_v_std_lambda_TUBE.set_xlabel(f'Zenith Angle {deg}')
# ax_zenith_v_std_lambda_TUBE.set_ylabel(f'Stadard Deviation of Wavelength\nAbsorbed by Tube {mum}')
# ax_zenith_v_std_lambda_TUBE.legend()
# ax_zenith_v_std_lambda_TUBE.grid(visible=True, alpha=.7)

ax_airmass_v_std_lambda_TUBE.set_title(f'Air Length vs. Stadard Deviation\nof Wavelength Absorbed by Tube')
ax_airmass_v_std_lambda_TUBE.set_xlabel(f'Air Length [m]')
ax_airmass_v_std_lambda_TUBE.set_ylabel(f'Stadard Deviation of Wavelength\nAbsorbed by Tube {mum}')
ax_airmass_v_std_lambda_TUBE.legend()
ax_airmass_v_std_lambda_TUBE.grid(visible=True, alpha=.7)
ax_airmass_v_std_lambda_TUBE.set_xscale('log')

# ax_scatters_v_std_lambda_TUBE.set_title(f'Mean Scattering Events vs. Standard\nDeviation of Wavelgth Absorbed by Tube')
# ax_scatters_v_std_lambda_TUBE.set_xlabel(f'Mean Scattering Events [-]')
# ax_scatters_v_std_lambda_TUBE.set_ylabel(f'Standard Deviation of Wavelgth\nAbsorbed by Tube {mum}')
# ax_scatters_v_std_lambda_TUBE.legend()
# ax_scatters_v_std_lambda_TUBE.grid(visible=True, alpha=.7)

# ax_zenith_v_std_lambda_SKY.set_title(f'Zenith Angle vs. Stadard Deviation\nof Wavelength Exiting Sky')
# ax_zenith_v_std_lambda_SKY.set_xlabel(f'Zenith Angle {deg}')
# ax_zenith_v_std_lambda_SKY.set_ylabel(f'Stadard Deviation of\nWavelength Exiting Sky {mum}')
# ax_zenith_v_std_lambda_SKY.legend()
# ax_zenith_v_std_lambda_SKY.grid(visible=True, alpha=.7)

ax_airmass_v_std_lambda_SKY.set_title(f'Air Length vs. Stadard Deviation\nof Wavelength Exiting Sky')
ax_airmass_v_std_lambda_SKY.set_xlabel(f'Air Length [m]')
ax_airmass_v_std_lambda_SKY.set_ylabel(f'Stadard Deviation of\nWavelength Exiting Sky {mum}')
ax_airmass_v_std_lambda_SKY.legend()
ax_airmass_v_std_lambda_SKY.grid(visible=True, alpha=.7)
ax_airmass_v_std_lambda_SKY.set_xscale('log')

# ax_scatters_v_std_lambda_SKY.set_title(f'Mean Scattering Events vs. Standard\nDeviation of Wavelgth Exiting Sky')
# ax_scatters_v_std_lambda_SKY.set_xlabel(f'Mean Scattering Events [-]')
# ax_scatters_v_std_lambda_SKY.set_ylabel(f'Standard Deviation of\nWavelgth Exiting Sky {mum}')
# ax_scatters_v_std_lambda_SKY.legend()
# ax_scatters_v_std_lambda_SKY.grid(visible=True, alpha=.7)

# ax.legend()
# 


# ax.scatter(df['zeniths'], df['arr_TIME'])
# ax.set_title(f'{name} - {rays[:2]}k: time')



plt.show()