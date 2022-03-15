import numpy as np
import random
from math import *


def schedule(T_pre):  # cooling schedule
    return T_pre / (1+0.001*T_pre)


def accept(delta, t): return exp(-delta/t)  # Acceptance probability function


def n_iter(t0, t): return 1.5*(t0-t) + 10  # n_iter function


def change(A_list, constraint, step, steps):
    # sample elements to change the Height
    if step/steps < 0.2:
        num = 6
    elif step/steps < 0.5:
        num = 4
    elif step/steps < 0.85:
        num = 2
    else:
        num = 1
    ele = random.sample([_ for _ in range(len(A_list))], num)
    idx = [i for i in range(len(A_list)) if i not in ele]
    upper = constraint - sum(A_list[idx])
    A_copy = A_list.copy()
    while True:
        A_copy[ele] = np.random.uniform(0.001, upper, num)
        if A_copy.sum() <= constraint:
            break
    A_copy = np.round(A_copy, 3)
    return A_copy
