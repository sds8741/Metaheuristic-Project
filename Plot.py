import matplotlib.pyplot as plt
import numpy as np


def history(steps, record, best_record):
    x = np.arange(1, steps+1)
    plt.plot(x, record, 'b--', label='Displacement record')
    plt.plot(x, best_record, 'r-', label='Converge history')
    plt.xlabel('Steps')
    plt.ylabel('Displacement (m)')
    plt.title("Converge History")
    plt.grid(True)
    plt.legend()
    # plt.show()


def structure(coord, n, L, H):
    y = np.zeros(n)
    coord = coord.ravel()
    mean = np.mean(H)
    max_ = max(H)
    for i in range(n-1):
        xx = coord[i:i+2]
        yy = y[i:i+2] + H[i]
        plt.plot(xx, yy, 'orange')
        plt.fill_between(xx, yy, color='orange')
    plt.xlim(0, L)
    plt.ylim(-max_/2, max_ + mean/1.5)
    plt.xlabel('X direction (m)')
    plt.ylabel('Height (m)')
    plt.title('Beam Shape')
    # plt.show()


def plot(steps, record, best_record, coord, n, L, H):
    plt.figure(figsize=(14, 5))
    plt.subplot(1, 2, 1)
    history(steps, record, best_record)
    plt.subplot(1, 2, 2)
    structure(coord, n, L, H)
    plt.show()
