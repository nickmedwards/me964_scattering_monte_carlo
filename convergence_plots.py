import ast
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['font.size'] = 22

converge_n_rays = [100, 1000, 10000, 100000]

f = open('.\convergence.txt')
lns = {}
for l in f.readlines():
    k, arr = l[:-1].split(' : ')
    lns[k] = np.array(ast.literal_eval(arr))

f.close()

fig, ax = plt.subplots(3,1, sharex=True)
ax[0].errorbar(converge_n_rays, 100 * lns['mean_N_ABS_TUBE'], 100 * lns['std_N_ABS_TUBE'], label='Tube')
ax[0].errorbar(converge_n_rays, 100 * lns['mean_N_ABS_SKY'], 100 * lns['std_N_ABS_SKY'], label='Ground')
ax[0].errorbar(converge_n_rays, 100 * lns['mean_N_ABS_GROUND'], 100 * lns['std_N_ABS_GROUND'], label='Sky')
ax[0].errorbar(converge_n_rays, 100 * lns['mean_N_ABS_TROUGH'], 100 * lns['std_N_ABS_TROUGH'], label='Trough')
# ax[0].set_xscale('log')
ax[0].legend()
# ax[0].set_title('Percent of rays absorbed by each geometry')
# ax[0].set_xlabel('Number of rays in simualtion [-]')
ax[0].set_ylabel('Hit Percent [%]')

# fig, ax = plt.subplots(1,1)
ax[1].errorbar(converge_n_rays, 100 * lns['mean_E_TUBE'], 100 * lns['std_E_TUBE'], label='Tube')
ax[1].errorbar(converge_n_rays, 100 * lns['mean_E_SKY'], 100 * lns['std_E_SKY'], label='Ground')
ax[1].errorbar(converge_n_rays, 100 * lns['mean_E_GROUND'], 100 * lns['std_E_GROUND'], label='Sky')
ax[1].errorbar(converge_n_rays, 100 * lns['mean_E_TROUGH'], 100 * lns['std_E_TROUGH'], label='Trough')
# ax[1].set_xscale('log')
# ax[1].legend()
# ax[1].set_title('Percent of total energy absorbed by each geometry')
# ax[1].set_xlabel('Number of rays in simualtion [-]')
ax[1].set_ylabel('Percent total energy [%]')

# fig, ax = plt.subplots(1,1)
ax[2].errorbar(converge_n_rays, lns['mean_E_TOTAL'], lns['std_E_TOTAL'], label='Total')
ax[2].set_xscale('log')
# ax[2].set_title('Average energy of ray vs number of rays in simulation')
ax[2].set_xlabel('Number of rays in simualtion [-]')
ax[2].set_ylabel('Average energy\nof ray [J]')

plt.show()