import math
import numpy as np
from constants import *
import numpy.linalg as linalg


def T1(theta, fy):
    """
    transition matrix from B to N frame
    :param theta:
    :param fy:
    :return:
    """
    return np.array([[np.cos(theta), 0, np.sin(theta)],
                     [np.sin(fy) * np.sin(theta), np.cos(fy), -np.sin(fy) * np.cos(theta)],
                     [-np.cos(fy) * np.sin(theta), np.sin(fy), np.cos(fy) * np.cos(theta)]])


def inverse_T1(theta, fy):
    return linalg.inv(T1(theta, fy))


def T2(betta):
    """
    transition matrix from D to B Frame
    :param betta:
    :return:
    """
    return np.array([[np.cos(betta), -np.sin(betta), 0],
                     [np.sin(betta), np.cos(betta), 0],
                     [0, 0, 1]])


def inverse_T2(betta):
    """ Transition matrix from B to D frame"""
    return linalg.inv(T2(betta))


def angular_velocity_vector_B_frame(fy_dot, theta, theta_dot, gamma_dot):
    return np.array([[fy_dot * np.cos(theta)],
                     [theta_dot],
                     [fy_dot * np.sin(theta) + gamma_dot]])


def F_d(alpha, V_x, V_y, V_z, theta, fy, betta):
    """
    returns F total in vector form in d frame of reference
    :param alpha:
    :param V:
    :return:
    """
    cl = calc_Cl(alpha)
    cd = calc_Cd(alpha)
    F = coef * calc_v_rel_3D(V_x, V_y, V_z) * np.array([[cl * np.sin(alpha) - cd * np.cos(alpha)],
                                                        [0],
                                                        [cl * np.cos(alpha) + cd * np.sin(alpha)]])  # in d frame
    # adding gravity
    inv_mat = inverse_T1(theta, fy)
    gravity_N = np.array([[0], [0], [m * g]])
    gravity_B = np.matmul(inv_mat, gravity_N)
    inv_mat = inverse_T2(betta)
    gravity_D = np.matmul(inv_mat, gravity_B)
    F = F + gravity_D
    return F


def M_d(gamma_dot, omega_B, fy_dot, theta, alpha, alpha_dot, betta, V_x, V_y, V_z):
    """moment in d frame"""
    omega_D = np.matmul(inverse_T2(betta), omega_B)
    o1 = omega_D[1][0]
    m1 = 0.5 * A * d * Ro * calc_v_rel_3D(V_x, V_y, V_z) * (CR_gamma_tag * gamma_dot + CM_alpha_tag * o1)
    m2 = 0.5 * A * d * Ro * calc_v_rel_3D(V_x, V_y, V_z) * (CM_alpha * alpha + CM_alpha_tag * alpha_dot)
    m3 = CNr * (fy_dot * np.sin(theta) + gamma_dot)
    return np.array([[m1], [m2], [m3]])


def vec_from_D_to_B(vec, betta):
    return np.matmul(T2(betta), vec)


def vec_from_N_to_B(vec, theta, fy):
    return np.matmul(T1(theta, fy), vec)


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2, z):
    """ Returns the angle in radians between vectors 'v1' and 'v2'"""
    v1_size = np.linalg.norm(v1)
    v2_size = np.linalg.norm(v2)
    if z == 17.91988015645636:
        print(v1, v2)
        print(v1_size, v2_size)
        print(np.vdot(v1, v2) / (v2.size * v1.size))
    if np.vdot(v1, v2) / (v2.size * v1.size) > 1:
            return np.arccos(1)
    return np.arccos(np.vdot(v1, v2) / (v2.size * v1.size))


def simulation(x_0=0, y_0=0, z_0=0, V0_x=0, V0_y=0, V0_z=0, theta0=0, theta_dot0=0, fy0=0, fy_dot0=0,
               gamma0=0, gamma_dot0=0):
    """עד מתייייייייייי"""
    X = [x_0]
    Y = [y_0]
    Z = [z_0]
    V_x = [V0_x]  # N frame of reference
    V_y = [V0_y]
    V_z = [V0_z]

    theta = [theta0]
    theta_dot = [theta_dot0]
    gamma = [gamma0]
    gamma_dot = [gamma_dot0]
    fy = [fy0]
    fy_dot = [fy_dot0]

    V_N = np.array([[V_x[-1]],
                    [V_y[-1]],
                    [V_z[-1]]])
    V_b = vec_from_N_to_B(V_N, theta[-1], fy[-1])
    n1 = np.array([[0], [0], [1]])
    b1 = vec_from_N_to_B(n1, theta[-1], fy[-1])
    alpha0 = angle_between(V_N, b1, Z[-1]) - np.pi / 2
    alpha = [alpha0]

    # B frame of reference
    U_x = [V_b[0][0]]
    U_y = [V_b[1][0]]
    U_z = [V_b[2][0]]

    while Z[-1] > 0:
        # print(Z[-1])
        V_N = np.array([[V_x[-1]],
                        [V_y[-1]],
                        [V_z[-1]]])
        V_b = vec_from_N_to_B(V_N, theta[-1], fy[-1])
        betta = np.arctan(V_b[0][0] / V_b[1][0])

        n1 = np.array([[0], [0], [1]])
        b1 = vec_from_N_to_B(n1, theta[-1], fy[-1])
        alpha.append(angle_between(V_N, b1,Z[-1]) - np.pi / 2)

        Fd = F_d(alpha[-1], V_x[-1], V_y[-1], V_z[-1], theta[-1], fy[-1], betta)
        Fb = vec_from_D_to_B(Fd, betta)

        md = M_d(gamma_dot[-1], angular_velocity_vector_B_frame(fy_dot[-1], theta[-1], theta_dot[-1], gamma_dot[-1]),
                 fy_dot[-1],
                 theta[-1], alpha[-1], (alpha[-1] - alpha[-2]) / dt, betta, V_x[-1], V_y[-1], V_z[-1])
        Mb = vec_from_D_to_B(md, betta)

        ux_tag = (1 / m) * (Fb[0][0])
        # - V_b[2] * theta_dot[-1] + V_b[1] * fy_dot[-1] * np.sin(theta[-1])
        U_x.append(U_x[-1] + dt * ux_tag)

        uy_tag = (1 / m) * (Fb[1][0])
        #  + V_b[2] * fy_dot[-1] * np.cos(theta[-1]) - V_b[0] * fy_dot[-1] * np.sin(theta[-1])
        U_y.append(U_y[-1] + dt * uy_tag)

        uz_tag = (1 / m) * (Fb[2][0])
        # + V_b[0] * theta_dot[-1] - V_b[1] * fy_dot[-1] * np.cos(theta[-1])
        U_z.append(U_y[-1] + dt * uz_tag)

        Unew = np.array([[U_x[-1]],
                         [U_y[-1]],
                         [U_z[-1]]])
        Vnew = np.matmul(T1(theta[-1], fy[-1]), Unew)

        V_x.append(Vnew[0][0])
        V_y.append(Vnew[1][0])
        V_z.append(Vnew[2][0])

        X.append(Vnew[0][0] * dt + X[-1])
        Y.append(Vnew[1][0] * dt + Y[-1])
        Z.append(Vnew[2][0] * dt + Z[-1])

        fy_dot2 = (1 / I_X_Y) * (Mb[0][0] + I_X_Y * theta_dot[-1] * fy_dot[-1] * np.sin(theta[-1])
                                 - I_z * theta_dot[-1] * (fy_dot[-1] * np.sin(theta[-1]) + gamma_dot[-1])
                                 + I_X_Y * theta_dot[-1] * fy_dot[-1] * np.sin(theta[-1])) * np.cos(theta[-1])
        theta_dot2 = (1 / I_X_Y) * (
                Mb[1][0] + I_z * fy_dot[-1] * np.cos(theta[-1]) * (fy_dot[-1] * np.sin(theta[-1]) + gamma_dot[-1])
                - I_x * fy_dot[-1] * fy_dot[-1] * np.sin(theta[-1]) * np.cos(theta[-1]))
        gamma_dot2 = (1 / I_z) * (
                Mb[2][0] - I_z * fy_dot2 * np.sin(theta[-1]) + theta_dot[-1] * fy_dot[-1] * np.cos(theta[-1]))

        gamma_dot.append(gamma_dot[-1] + dt * gamma_dot2)
        gamma.append(gamma[-1] + dt * gamma_dot[-1])

        theta_dot.append(theta_dot[-1] + dt * theta_dot2)
        theta.append(theta[-1] + dt * theta_dot[-1])

        fy_dot.append(fy[-1] + dt * fy_dot2)
        fy.append(fy[-1] + dt * fy_dot[-1])

    return X, Y, Z, V_x, V_y, V_z, theta, theta_dot, fy, fy_dot, gamma, gamma_dot


if __name__ == '__main__':
    X, Y, Z, V_x, V_y, V_z, theta, theta_dot, fy, fy_dot, gamma, gamma_dot = simulation(
        0., 0., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0.
    )
   #  print(Z)
    draw_3D_graph(X, 'x', Y, 'y', Z, 'z')
    # draw_xyt(x, 'x', z, 'z', t, 't')
    # draw_y_as_x(t, 't', ux, 'u_x')
    # draw_y_as_x(t, 't', uz, 'u_z')
    # draw_y_as_x(t, 't', np.array(theta) * 180 / np.pi, 'theta')
