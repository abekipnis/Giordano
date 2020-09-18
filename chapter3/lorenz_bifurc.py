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

    d = pd.read_csv(args.files[0])
    plt.scatter(d.r, d.z, s=0.02)
    plt.xlabel("r parameter")
    plt.ylabel("z values of local maxima\nafter transients have died away")
    plt.show()

    # d = [pd.read_csv(f) for f in args.files]
    # plt.title("Bifurcation diagram for driven pendulum")
    #
    # scatter = [[],[]] #driving force vs maxima
    # for n in range(len(d)):
    #     r = re.search(r'r-(.*?).dat', args.files[n]).group(1)
    #     for dp in range(len(d[n])):
    #         scatter[0].append(r)
    #         scatter[1].append(d[n].iloc[dp].z)
    #
    # plt.scatter(scatter[0], abs(scatter[1]), s=0.05, c='black')
    #
    #
    # ax = plt.gca()
    # start, end = ax.get_xlim()
    #
    # ax.xaxis.set_ticks(np.arange(start, end, 0.15))
    # plt.xlabel("Driving force")
    # plt.ylabel("Fixed point")
    # plt.show()
