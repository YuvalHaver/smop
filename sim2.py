import math
import numpy as np
import freesbe_graphs as fg

ro = 1.2041
cl0 = 0.15
cl_alpha = 1.4
cd0 = 0.08
cd_alpha = 2.72
alpha0 = -4

R = 0.15
g = 9.81
m = 0.175
A = math.pi * R ** 2
d = 2 * R

coef = (ro * A) / (2 * m)
dt = 10 ** -4

I_x = 0.5 * m * R ** 2
CRr = 0.0017
CRp = -5.5 * 10 ** (-3)
CM0 = -0.08
CM_alpha = 0.43
CMq = -5 * 10 ** (-3)
CNr = -7.18 * 10 ** (-6)


def calculate_Cl(alpha):
    cl = cl0 + cl_alpha * alpha * 2 * math.pi / 360
    return cl


def calculete_Fl(cl, v):
    F = 0.5 * ro * A * (v ** 2) * cl
    return F


def calculate_Cd(alpha):
    cd = cd0 + cd_alpha * ((alpha - alpha0) * 2 * math.pi / 360) ** 2
    return cd


def calculate_Fd(cd, v):
    F = 0.5 * ro * A * (v ** 2) * cd


def calculate_v_rel(ux, uz):
    return np.sqrt(math.pow(ux, 2) + math.pow(uz, 2))


def one_step(u_x_prev, u_z_prev, x_arr, z_arr, theta):
    betta = np.arctan(float(u_z_prev[-1]) / float(u_x_prev[-1]))
    alpha = betta / theta[-1]
    cl = calculate_Cl(alpha)
    cd = calculate_Cd(alpha)
    ux_new = u_x_prev[-1] + dt * coef * math.pow(calculate_v_rel(u_x_prev[-1], u_z_prev[-1]), 2) * (
                -cl * np.sin(betta) - cd * np.cos(betta))
    uz_new = u_z_prev[-1] + dt * (coef * math.pow(calculate_v_rel(u_x_prev[-1], u_z_prev[-1]), 2) * (
                cl * np.cos(betta) - cd * np.sin(betta)) - g)
    x_arr.append(x_arr[-1] + ux_new * dt)
    z_arr.append(z_arr[-1] + uz_new * dt)

    betta = np.arctan(float(u_z_prev[-1]) / float(u_x_prev[-1]))

    theta.append(((CM0 + CM_alpha * (theta[-1] - betta) + CMq * ((theta[-1] - theta[-2]) / dt)) * 0.5 * ro * (
                calculate_v_rel(u_x_prev[-1], u_z_prev[-1]) ** 2) * A * (dt ** 2) * d / I_x) + 2 * theta[-1] - theta[
                     -2])
    return ux_new, uz_new


def basic_simulation(u0_x, u0_z, theta, x_0, z_0):
    theta = [theta, theta]
    x_array = [x_0]
    z_array = [z_0]
    ux_arr = [u0_x]
    uz_arr = [u0_z]
    u_x_curr = u0_x
    u_z_curr = u0_z
    t = [0]

    while (z_array[-1] > 0):
        u_x_curr, u_z_curr = one_step(ux_arr, uz_arr, x_array, z_array, theta)
        ux_arr.append(u_x_curr)
        uz_arr.append(u_z_curr)
        t.append(t[-1] + dt)
    theta.remove(theta[-1])
    return x_array, z_array, t, ux_arr, uz_arr, theta


if __name__ == '__main__':
    x, z, t, ux, uz, theta = basic_simulation(10, 0, 45, 0, 200)
    fg.draw_y_as_x(t, 't', theta, 'theta')
    fg.draw_xyt(x, 'x', z, 'z', t, 't')
    fg.draw_y_as_x(t, 't', ux, 'u_x')
    fg.draw_y_as_x(t, 't', uz, 'u_z')
    betta_arr = [np.arctan(uz[i] / ux[i]) * (180 / np.pi) for i in range(len(ux))]
    fg.draw_y_as_x(t, 't', betta_arr, 'betta')
