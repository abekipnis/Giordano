import pandas as pd
import matplotlib.pyplot as plt
import argparse
import re
import numpy as np

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--files', type=str,  nargs='+', help='file(s) to read')
    # parser.add_argument('--driving_freq', type=float)
    args = parser.parse_args()

    d = [pd.read_csv(f) for f in args.files]
    plt.title("Bifurcation diagram for driven pendulum")

    scatter = []
    for n in range(len(d)):
        f_d = re.search(r'f_d-(.*?).dat', args.files[n]).group(1)
        f = [[float(f_d)]*5, d[n].theta.iloc[-5:]]
        scatter.append(np.array(f).T)

    scatter = np.array(scatter).T
    plt.scatter(scatter[0].flatten(), abs(scatter[1].flatten()), s=0.05, c='black')


    ax = plt.gca()
    start, end = ax.get_xlim()

    ax.xaxis.set_ticks(np.arange(start, end, 0.15))
    plt.xlabel("Driving force")
    plt.ylabel("Fixed point")
    # plt.xticks([1.35, 1.4, 1.45, 1.5],[1.35, 1.4, 1.45, 1.5])
    plt.show()
