import math
import numpy as np
from constants import *


def one_step(u_x_prev, u_y_prev, u_z_prev, x_arr, y_arr, z_arr, theta, theta_1, theta_2, phi, phi_1, phi_2, gamma, gamma_1, gamma_2):
    betta = np.arctan(float(u_z_prev[-1]) / float(u_x_prev[-1]))
    v_rel=calc_v_rel_3D(u_x_prev[-1], u_y_prev[-1], u_z_prev[-1])
    alpha = (theta[-1] - betta)
    cl = calc_Cl(alpha)
    cd = calc_Cd(alpha)
    u_x_prev.append(u_x_prev[-1] + dt * coef * math.pow(v_rel, 2) * (
            -cl * np.sin(betta) - cd * np.cos(betta)))
    u_y_prev.append(u_y_prev[-1] + dt * coef * math.pow(v_rel, 2) * (
            -cl * np.sin(betta) - cd * np.cos(betta)))
    u_z_prev.append(u_z_prev[-1] + dt * (coef * math.pow(v_rel, 2) * (
            cl * np.cos(betta) - cd * np.sin(betta)) - g))
    x_arr.append(x_arr[-1] + u_x_prev[-1] * dt)
    y_arr.append(y_arr[-1] + u_y_prev[-1] * dt)
    z_arr.append(z_arr[-1] + u_z_prev[-1] * dt)


def basic_simulation(u0_x, u0_y, u0_z, theta, x_0, y_0, z_0):
    gamma = [0]#z axis
    gamma_1 = [0]
    gamma_2 = [0]
    phi = [0]#y axis
    phi_1 = [0]
    phi_2 = [0]
    theta = [theta]#x axis
    theta_1 = [0]
    theta_2 = [0]
    x_array = [x_0]
    y_array=[y_0]
    z_array = [z_0]
    ux_arr = [u0_x]
    uy_arr = [u0_y]
    uz_arr = [u0_z]
    t = [0]

    while (z_array[-1] > 0):
        one_step(ux_arr, uy_arr, uz_arr, x_array, y_array, z_array, theta, theta_1, theta_2, phi, phi_1, phi_2, gamma, gamma_1, gamma_2)
        t.append(t[-1] + dt)
    return ux_arr, uy_arr, uz_arr, theta, x_array, y_array, z_array, t


if __name__ == '__main__':
    ux, uy, uz, theta, x, y, z, t = basic_simulation(14, 0, 0,0 * np.pi / 180, 0, 0, 2)
    # draw_xyt(x, 'x', z, 'z', t, 't')
    draw_y_as_x(x, 'X', z, 'Z')
    # draw_y_as_x(t, 't', ux, 'u_x')
    # draw_y_as_x(t, 't', uz, 'u_z')
    draw_y_as_x(t, 't', np.array(theta) * 180 / np.pi, 'theta')
    # comment
