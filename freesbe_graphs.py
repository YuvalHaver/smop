import numpy as np
import matplotlib.pyplot as plt

def draw_y_as_x (x,x_name, y, y_name):
    """
    draw graph y(x)
    :param x: horizontal values array
    :param x_name: name of x axis
    :param y: vertical values array
    :param y_name: name of y axis
    :return: None
    """
    plt.plot(x, y)
    plt.title(f'{y_name}({x_name})')
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.grid()
    plt.show()

def draw_xyt(x,x_name, y, y_name, t, t_name):
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
    draw_y_as_x(x, x_name, t, t_name)
    draw_y_as_x(y, y_name, t, t_name)
    draw_y_as_x(y, y_name, x, x_name)
