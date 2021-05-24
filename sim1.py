import math
import numpy as np
import freesbe_graphs as fg

ro = 1.2041
cl0 = 0.15
cl_alpha = 2.8
cd0 = 0.1
cd_alpha = 2.5
alpha0=-4

g=9.81
m=0.175
A = math.pi * 0.15 ** 2

coef=(ro*A)/(2*m)
dt=10**-4


def calculate_Cl(alpha):
    cl = cl0 + cl_alpha * alpha * 2 * math.pi / 360
    return cl


def calculete_Fl(cl, v):
    F = 0.5 * ro * A * (v ** 2) * cl
    return F

def calculate_Cd(alpha):
    cd=cd0+cd_alpha*((alpha-alpha0)*2 * math.pi / 360) **2
    return cd

def calculate_Fd(cd,v):
    F = 0.5 * ro * A * (v ** 2) * cd

def calculate_v_rel(ux,uz, vx,vz):

    return np.sqrt(math.pow(ux,2)+math.pow(uz,2))

def one_step(u_x_prev, u_z_prev,cl,cd,vx,vz,x_arr,z_arr):

    betta=np.arctan(float(u_z_prev[-1])/float(u_x_prev[-1]))
    ux_new=u_x_prev[-1]+dt*coef*math.pow(calculate_v_rel(u_x_prev[-1],u_z_prev[-1],vx,vz),2)*(-cl*np.sin(betta)-cd*np.cos(betta))
    uz_new=u_z_prev[-1]+dt*(coef*math.pow(calculate_v_rel(u_x_prev[-1],u_z_prev[-1],vx,vz),2)*(cl*np.cos(betta)-cd*np.sin(betta))-g)
    x_arr.append(x_arr[-1]+ux_new*dt)
    z_arr.append(z_arr[-1]+uz_new*dt)
    return ux_new,uz_new


def basic_simulation(u0_x,u0_z,v_x_air,v_z_air,alpha,x_0,z_0):

    x_array=[x_0]
    z_array=[z_0]
    ux_arr=[u0_x]
    uz_arr=[u0_z]
    u_x_curr=u0_x
    u_z_curr=u0_z
    t = [0]

    cl=calculate_Cl(alpha)
    cd=calculate_Cd(alpha)
    while (z_array[-1]>0):
        u_x_curr,u_z_curr=one_step(ux_arr,uz_arr,cl,cd,v_x_air,v_z_air, x_array,z_array)
        ux_arr.append(u_x_curr)
        uz_arr.append(u_z_curr)
        t.append(t[-1] + dt)

    return x_array, z_array, t, ux_arr, uz_arr

if __name__ == '__main__':
    x,z, t,ux,uz = basic_simulation(14,0,0,0,10,0,1)
    fg.draw_xyt(x, 'x', z, 'z', t, 't')
    fg.draw_y_as_x( t, 't', ux, 'u_x')
    fg.draw_y_as_x(  t, 't', uz, 'u_z')
    betta_arr=[np.arctan(uz[i]/ux[i])*(180/np.pi) for i in range(len(ux))]
    fg.draw_y_as_x(t,'t', betta_arr, 'betta')



