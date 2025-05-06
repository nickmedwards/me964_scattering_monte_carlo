from monte_carlo import *

test_trough = Trough(W_T, L_T, M_PARABOLA[1], RHO_PARABOLA[1])
tube = Tube(L_T, 1, ALPHA_TUBE[1])
sky = Sky()
sun = Sun(np.pi / 8, test_trough)

######################
# tube hit test rays #
######################

# sky_ray = Ray(
#     .5,
#     np.array([0, 0., -5]),
#     sph2rec(1, 0.00000001, -np.pi / 4)
# ) # ccw 
# hit_trace(sky_ray, sky, test_trough, tube)

# # cw
# test_ray = Ray(
#     5,
#     np.array([5, 0., 5]),
#     sph2rec(1, 0.00000001, -2 * np.pi / 3)
# )
# hit_trace(test_ray, sky, test_trough, tube)

# # cw
# test_ray = Ray(
#     5,
#     np.array([5, 0., 5]),
#     np.array([-.5, -8.66025404e-09, -8.66025404e-01])
# )
# hit_trace(test_ray, sky, test_trough, tube)

# # perfect hit
# test_ray = Ray(
#     5,
#     np.array([5, 0., 5]),
#     sph2rec(1, 0.00000001, -3 * np.pi / 4)
# )
# hit_trace(test_ray, sky, test_trough, tube)