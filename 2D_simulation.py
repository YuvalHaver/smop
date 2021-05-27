import math
import numpy as np
import constants as cst


def calc_next_step(V_x_prev, V_z_prev, theta, X, Z):
    betta = np.arctan(float(V_z_prev) / float(V_x_prev))
    cl = cst.calc_Cl(theta - betta)
    cd = cst.calc_Cd(theta - betta)
    V_x_new = V_x_prev + cst.dt * cst.coef * math.pow(cst.calc_v_rel(V_x_prev, V_z_prev), 2) * (-cl * np.sin(betta) - cd * np.cos(betta))
    V_z_new = V_z_prev + cst.dt * (cst.coef * math.pow(cst.calc_v_rel(V_x_prev, V_z_prev), 2) * (cl * np.cos(betta) - cd * np.sin(betta)) - cst.g)
    X.append(X[-1] + V_x_new * cst.dt)
    Z.append(Z[-1] + V_z_new * cst.dt)
    return V_x_new, V_z_new


def basic_simulation(V0_x, V0_z, theta, X_0, Z_0):
    X = [X_0]
    Z = [Z_0]
    V_x = [V0_x]
    V_z = [V0_z]
    V_x_curr = V0_x
    V_z_curr = V0_z
    t = [0]


    while Z[-1] > 0:
        V_x_curr, V_z_curr = calc_next_step(V_x[-1], V_z[-1], theta, X, Z)
        V_x.append(V_x_curr)
        V_z.append(V_z_curr)
        t.append(t[-1] + cst.dt)

    return X, Z, t, V_x, V_z


if __name__ == '__main__':
    x, z, t, ux, uz = basic_simulation(10, 0, 7 * np.pi/180, 0, 2)
    cst.draw_xyt(x, 'x', z, 'z', t, 't')
    cst.draw_y_as_x(t, 't', ux, 'u_x')
    cst.draw_y_as_x(t, 't', uz, 'u_z')
    betta_arr = [np.arctan(uz[i] / ux[i]) * (180 / np.pi) for i in range(len(ux))]
    cst.draw_y_as_x(t, 't', betta_arr, 'betta')
