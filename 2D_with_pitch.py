import math
import numpy as np
from constants import *


def one_step(u_x_prev, u_z_prev, x_arr, z_arr, theta, theta_1, theta_2, phi, phi_1, phi_2, gamma, gamma_1, gamma_2):
    betta = np.arctan(float(u_z_prev[-1]) / float(u_x_prev[-1]))
    alpha = (theta[-1] - betta)
    cl = calc_Cl(alpha)
    cd = calc_Cd(alpha)
    u_x_prev.append(u_x_prev[-1] + dt * coef * math.pow(calc_v_rel(u_x_prev[-1], u_z_prev[-1]), 2) * (
            -cl * np.sin(betta) - cd * np.cos(betta)))
    u_z_prev.append(u_z_prev[-1] + dt * (coef * math.pow(calc_v_rel(u_x_prev[-1], u_z_prev[-1]), 2) * (
            cl * np.cos(betta) - cd * np.sin(betta)) - g))
    x_arr.append(x_arr[-1] + u_x_prev[-1] * dt)
    z_arr.append(z_arr[-1] + u_z_prev[-1] * dt)

    betta = np.arctan(float(u_z_prev[-1]) / float(u_x_prev[-1]))

    constants = 0.5 * Ro * (
            calc_v_rel(u_x_prev[-1], u_z_prev[-1]) ** 2) * A * (dt ** 2) * d

    pitch = (CM0 + CM_alpha * (theta[-1] - betta) + CMq * (theta_1[-1])) * constants

    theta_2.append((pitch + I_z * phi_1[-1] * np.cos(theta[-1]) * (phi_1[-1] * np.sin(theta[-1]) + gamma_1[-1]) - I_z *
                    phi_1[-1] ** 2 * np.cos(theta[-1]) * np.sin(theta[-1])) / I_x)

    theta_1.append(theta_1[-1] + theta_2[-1] * dt)

    theta.append((theta[-1] + dt * theta_1[-1]) % (np.pi * 2))

    roll = (CRr * gamma_1[-1] + CRp * phi_1[-1]) * constants

    phi_2.append((roll + I_x * theta_1[-1] * phi_1[-1] * np.sin(theta[-1]) - I_z * theta_1[-1] * (
                phi_1[-1] * np.sin(theta[-1]) + gamma_1[-1]) + I_x * theta_1[-1] * phi_1[-1]) * np.cos(theta[-1]) / I_x)

    phi_1.append(phi_1[-1] + phi_2[-1] * dt)

    phi.append((phi[-1] + dt * phi_1[-1]) % (np.pi * 2))

    spin_down = (CNr * gamma_1[-1]) * constants

    gamma_2.append((spin_down - I_z * phi_1[-1] * np.sin(theta[-1]) + theta_1[-1] * phi_1[-1] * np.cos(theta[-1])) / I_z)
    gamma_1.append(gamma_1[-1] + gamma_2[-1] * dt)

    gamma.append((gamma[-1] + dt * gamma_1[-1]) % (np.pi * 2))

    # theta.append(((CM0 + CM_alpha * (theta[-1] - betta) + CMq * ((theta[-1] - theta[-2]) / dt)) * 0.5 * Ro * (
    #         calc_v_rel(u_x_prev[-1], u_z_prev[-1]) ** 2) * A * (dt ** 2) * d / I_x) + 2 * theta[-1] - theta[
    #                  -2])
    # omega.append(omega[-1] + dt / I_z * CNr * omega[-1] * 0.5 * Ro * A * d * (
    #         calc_v_rel(u_x_prev[-1], u_z_prev[-1]) ** 2))


def basic_simulation(u0_x, u0_z, theta, x_0, z_0):
    gamma = [0]
    gamma_1 = [0]
    gamma_2 = [0]
    phi = [0]
    phi_1 = [0]
    phi_2 = [0]
    theta = [theta]
    theta_1 = [0]
    theta_2 = [0]
    x_array = [x_0]
    z_array = [z_0]
    ux_arr = [u0_x]
    uz_arr = [u0_z]
    t = [0]

    while (z_array[-1] > 0):
        one_step(ux_arr, uz_arr, x_array, z_array, theta, theta_1, theta_2, phi, phi_1, phi_2, gamma, gamma_1, gamma_2)
        t.append(t[-1] + dt)
    return  ux_arr, uz_arr, theta, x_array, z_array, t

if __name__ == '__main__':
    ux, uz, theta, x, z , t = basic_simulation(14, 0, 0 * np.pi/180, 0, 2)
    draw_xyt(x, 'x', z, 'z', t, 't')
    draw_y_as_x(t, 't', np.array(theta) * 180 / np.pi, 'theta')
    # comment
