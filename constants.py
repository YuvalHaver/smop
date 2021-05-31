import numpy as np
import matplotlib.pyplot as plt
import math

"constants"

Cl_0 = 0.15
Cl_alpha = 1.41
cd_0 = 0.08
cd_alpha = 2.7
alpha_0 = -4 * (np.pi / 180)  # degrees

g = 9.81
Ro = 1.2041
m = 0.175
R = 0.1
d = 2*R
A = np.pi * R ** 2

coef = (Ro * A) / (2 * m)
dt = 10 ** -4

I_x = 0.25 * m * R ** 2
I_z = 0.5 * m * R ** 2
I_X_Y= -m*R**2

CRr = 0.0017
CRp = -5.5 * 10 ** (-3)
CM0 = -0.08
CM_alpha = 0.43
CMq = -5 * 10 ** (-3)
CNr = -7.18 * 10 ** (-6)


def calc_Cl(alpha):
    return Cl_0 + Cl_alpha * alpha


def calc_Cd(alpha):
    return cd_0 + cd_alpha * (alpha - alpha_0) ** 2


def calc_v_rel(V_x, V_z):
    return np.sqrt(math.pow(V_x, 2) + math.pow(V_z, 2))

def calc_v_rel_3D(V_x,V_y, V_z):
    return np.sqrt(math.pow(V_x, 2)+ math.pow(V_y, 2) + math.pow(V_z, 2))

def draw_y_as_x(x, x_name, y, y_name):
    """
    draw graph y(x)
    :param x: horizontal values array
    :param x_name: name of x axis
    :param y: vertical values array
    :param y_name: name of y axis
    :return: None
    """
    plt.scatter(x, y, s=0.1)
    plt.title(f'{y_name}({x_name})')
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.grid()
    plt.show()


def draw_xyt(x, x_name, y, y_name, t, t_name):
    """
    draws 3 graphs y(t), x(t), y(x)
    :param x: horizontal values array
    :param x_name: name of x axis
    :param y: vertical values array
    :param y_name: name of y axis
    :param t: t values
    :param t_name: t name
    :return: None
    """
    draw_y_as_x(t, t_name, x, x_name)
    draw_y_as_x(t, t_name, y, y_name)
    draw_y_as_x(x, x_name, y, y_name)


def draw_3D_graph(x, x_name, y, y_name, z, z_name):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(x, y, z)
    ax.set_xlabel(x_name)
    ax.set_ylabel(y_name)
    ax.set_zlabel(z_name)
    fig.show()
