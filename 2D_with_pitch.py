import math
import numpy as np
import freesbe_graphs as fg

ro = 1.2041
cl0 = 0.15
cl_alpha = 1.4
cd0 = 0.08
cd_alpha = 2.72
alpha0 = -4*np.pi/180

R = 0.15
g = 9.81
m = 0.175
A = math.pi * R ** 2
d = 2 * R

coef = (ro * A) / (2 * m)
dt = 10 ** -2

I_x = 0.25 * m * R ** 2
I_z = 0.5 * m * R ** 2

CRr = 0.0017
CRp = -5.5 * 10 ** (-3)
CM0 = -0.08
CM_alpha = 0.43
CMq = -5 * 10 ** (-3)
CNr = -7.18 * 10 ** (-6)


def calculate_Cl(alpha):
    cl = cl0 + cl_alpha * alpha
    return cl

def calculate_Cd(alpha):
    cd = cd0 + cd_alpha * ((alpha - alpha0)** 2)
    return cd

def calculate_v_rel(ux, uz):
    return np.sqrt(math.pow(ux, 2) + math.pow(uz, 2))


def one_step(u_x_prev, u_z_prev, x_arr, z_arr, theta, omega):
    betta = np.arctan(float(u_z_prev[-1]) / float(u_x_prev[-1]))
    alpha = theta[-1]-betta
    cl = calculate_Cl(alpha)
    cd = calculate_Cd(alpha)
    u_x_prev.append( u_x_prev[-1] + dt * coef * math.pow(calculate_v_rel(u_x_prev[-1], u_z_prev[-1]), 2) * (
            -cl * np.sin(betta) - cd * np.cos(betta)))
    u_z_prev.append(u_z_prev[-1] + dt * (coef * math.pow(calculate_v_rel(u_x_prev[-1], u_z_prev[-1]), 2) * (
            cl * np.cos(betta) - cd * np.sin(betta)) - g))
    x_arr.append(x_arr[-1] + u_x_prev[-1] * dt)
    z_arr.append(z_arr[-1] + u_z_prev[-1] * dt)

    betta = np.arctan(float(u_z_prev[-1]) / float(u_x_prev[-1]))

    theta.append(((CM0 + CM_alpha * (theta[-1] - betta) + CMq * ((theta[-1] - theta[-2]) / dt)) * 0.5 * ro * (
            calculate_v_rel(u_x_prev[-1], u_z_prev[-1]) ** 2) * A * (dt ** 2) * d / I_x) + 2 * theta[-1] - theta[
                     -2])

    # omega.append(omega[-1] + dt / I_z * CNr * omega[-1] * 0.5 * ro * A * d * (
    #         calculate_v_rel(u_x_prev[-1], u_z_prev[-1]) ** 2))


def basic_simulation(u0_x, u0_z, theta, omega, x_0, z_0):
    omega= [omega]
    theta = [theta, theta]
    x_array = [x_0]
    z_array = [z_0]
    ux_arr = [u0_x]
    uz_arr = [u0_z]
    t = [0]

    while (z_array[-1] > 0):
        one_step(ux_arr, uz_arr, x_array, z_array, theta, omega)
        t.append(t[-1] + dt)
    theta.remove(theta[-1])
    return x_array, z_array, t, ux_arr, uz_arr, theta, omega


if __name__ == '__main__':
    x, z, t, ux, uz, theta, omega = basic_simulation(14, 0, 14 / 180*np.pi, 3,  0, 1.7)
    fg.draw_xyt(x, 'x', z, 'z', t, 't')
    fg.draw_y_as_x(t, 't', ux, 'u_x')
    fg.draw_y_as_x(t, 't', uz, 'u_z')
    # comment
