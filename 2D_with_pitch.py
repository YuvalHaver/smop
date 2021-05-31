import math
import numpy as np
from constants import *

def one_step(u_x_prev, u_z_prev, x_arr, z_arr, theta, theta_1, theta_2, omega):
    betta = np.arctan(float(u_z_prev[-1]) / float(u_x_prev[-1]))
    alpha = (theta[-1]-betta)
    cl = calc_Cl(alpha)
    cd = calc_Cd(alpha)
    u_x_prev.append( u_x_prev[-1] + dt * coef * math.pow(calc_v_rel(u_x_prev[-1], u_z_prev[-1]), 2) * (
            -cl * np.sin(betta) - cd * np.cos(betta)))
    u_z_prev.append(u_z_prev[-1] + dt * (coef * math.pow(calc_v_rel(u_x_prev[-1], u_z_prev[-1]), 2) * (
            cl * np.cos(betta) - cd * np.sin(betta)) - g))
    x_arr.append(x_arr[-1] + u_x_prev[-1] * dt)
    z_arr.append(z_arr[-1] + u_z_prev[-1] * dt)

    betta = np.arctan(float(u_z_prev[-1]) / float(u_x_prev[-1]))
    theta_2.append((CM0 + CM_alpha * (theta[-1] - betta) + CMq * ((theta[-1] - theta[-2]) / dt)) * 0.5 * Ro * (
            calc_v_rel(u_x_prev[-1], u_z_prev[-1]) ** 2) * A * (dt ** 2) * d / I_x)
    theta_1.append(theta_1[-1]+theta_2[-1]*dt)
    theta.append((theta[-1]+dt*theta_1[-1])%(np.pi*2))
    # theta.append(((CM0 + CM_alpha * (theta[-1] - betta) + CMq * ((theta[-1] - theta[-2]) / dt)) * 0.5 * Ro * (
    #         calc_v_rel(u_x_prev[-1], u_z_prev[-1]) ** 2) * A * (dt ** 2) * d / I_x) + 2 * theta[-1] - theta[
    #                  -2])
    omega.append(omega[-1] + dt / I_z * CNr * omega[-1] * 0.5 * Ro * A * d * (
            calc_v_rel(u_x_prev[-1], u_z_prev[-1]) ** 2))
    print(omega[-1])

def basic_simulation(u0_x, u0_z, theta, omega, x_0, z_0):
    omega= [omega]
    theta = [theta, theta]
    theta_1=[0]
    theta_2=[0]
    x_array = [x_0]
    z_array = [z_0]
    ux_arr = [u0_x]
    uz_arr = [u0_z]
    t = [0]

    while (z_array[-1] > 0):
        one_step(ux_arr, uz_arr, x_array, z_array, theta,theta_1, theta_2, omega)
        t.append(t[-1] + dt)
    theta.remove(theta[-1])
    return x_array, z_array, t, ux_arr, uz_arr, theta, omega


if __name__ == '__main__':
    x, z, t, ux, uz, theta, omega = basic_simulation(14, 0, 10*np.pi/180, 3,  0, 1.7)
    draw_xyt(x, 'x', z, 'z', t, 't')
    # draw_y_as_x(t, 't', ux, 'u_x')
    # draw_y_as_x(t, 't', uz, 'u_z')
    draw_y_as_x(t, 't', np.array(theta)*180/np.pi, 'theta')
    # comment