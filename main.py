from math import inf
from random import random
import Solve
import SA
import Input
import numpy as np
import Plot


n = 15  # number of nodes (Odd number)
type_ = 1  # [ 1 for two fixed end beam ,others for cantilever]
force = 100  # concentrate load (kN)
load = 100  # weight load (kN/m)

total_length = 10  # Toatal length (m)
l = total_length / (n-1)  # element length
W = 0.5  # width(constant)
V = 1  # contraint volumn
constraint = (V/W)/l  # total height contraint (constraint total volumn)


# random initalize elements height list (m)
H = Input.random_H(n-1, constraint)

# SA algoritm parameter
SA_setting = {
    't0': 300,  # Intital Temperature
    'steps': 800,  # Iteration number
}

# Beam Parameter
config = {
    'NNOD': n,
    'NBC': n-1,
    'NMAT': 1,
    'NSEC': n-1,
    'NNE': 2,
    'IFORCE': 2,
    'COOR': np.linspace(0, total_length, n).reshape(1, -1),
    'NFIX': Input.NFIX(n, type_),
    'EXLD': Input.EXLD(n, force, type_),
    'IDBC': Input.IDBC(n-1),
    'PROP': np.array([[200e6], [0], [0], [0], [0]]),
    'SECT': Input.SECT(n-1, W, H),
    'FEF': Input.FEF(n-1, l, load),
}

# SA main code
dis_record, dis_record_best = [], []
best = inf
loc = n-3 if type_ == 1 else -2
for step in range(SA_setting['steps']):
    if step == 0:
        t = SA_setting['t0']
        Delta = Solve.Execute(**config)[loc].item()
    else:
        # cooling schedule
        t = SA.schedule(t)
    # n_iter function
    n_iter = SA.n_iter(SA_setting['t0'], t)
    for iter in range(round(n_iter)):
        H_next = SA.change(H, constraint, step, SA_setting['steps'])
        config['SECT'] = Input.SECT(n-1, W, H_next)
        Delta_next = Solve.Execute(**config)[loc].item()
        diff = abs(Delta_next) - abs(Delta)
        if diff < 0 or random() < SA.accept(diff*500, t):
            H = H_next
            break

    config['SECT'] = Input.SECT(n-1, W, H)

    Delta = Solve.Execute(**config)[loc].item()

    # save the best solution
    if abs(Delta) < best:
        H_best = H
        best = abs(Delta)

    dis_record.append(abs(Delta))
    dis_record_best.append(best)

    print(f'Step {step+1:3d}: Displacement {abs(Delta):.3f} m')

print('-------------------------')
print(f'Total volumn {(H_best*W*l).sum():.3f} m^3')
print('optimal H list', H_best)
print('optimal Displacement {:.3f} m'.format(best))
print(Solve.Execute(**config))
Plot.plot(SA_setting['steps'], dis_record, dis_record_best,
          config['COOR'], n, total_length, H_best)
