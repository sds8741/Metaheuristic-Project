import numpy as np
from math import floor
import random


def NFIX(n, type_):
    NFIX = np.zeros([2, n])
    if type_ == 1:
        NFIX[:, [0, -1]] = -1
    else:
        NFIX[:, 0] = -1
    return NFIX


def EXLD(n, force, type_):
    EXLD = np.zeros([2, n])
    if type_ == 1:
        node = floor(n/2)
        EXLD[0, int(node)] = -force
    else:
        EXLD[0, -1] = -force
    return EXLD


def IDBC(e):
    IDBC = np.zeros([5, e])
    IDBC[2, :] = 1
    IDBC[3, :] = np.arange(1, e+1)
    IDBC[0] = np.arange(1, e+1)
    IDBC[1] = IDBC[0] + 1
    return IDBC.astype('int')


def SECT(e, W, h_list):
    h_array = np.array(h_list)
    SECT = np.zeros([5, e])
    I_array = W * h_array**3 / 12
    SECT[1, :] = I_array
    return SECT


def FEF(e, l, load):
    FEF = np.zeros([4, e])
    FEF[[0, 2], :] = load*l / 2
    FEF[1, :] = load * l**2 / 12
    FEF[3, :] = -load * l**2 / 12
    return FEF


def random_H(e, constraint):
    upper = constraint / (e-1)
    while True:
        h_list = np.random.uniform(0.001, upper, e)
        if h_list.sum() <= constraint:
            break
    h_list = np.round(h_list, 3)
    return h_list
