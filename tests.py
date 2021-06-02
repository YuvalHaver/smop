import constants as cst
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from final_try import *
import numpy.linalg as linalg

if __name__ == '__main__':
    a1 = np.array([[1], [2], [3]])
    a3 = np.array([[1], [1], [1], [1]])
    a2 = np.mat([[1, 1, 1, 1], [1, 1, 1, 0], [1, 1, 0, 0], [1, 0, 0, 0]])
    #print(np.vdot(a3, a1))
    print(linalg.norm(a1))
    print(np.matmul(linalg.inv(a2), a2))
