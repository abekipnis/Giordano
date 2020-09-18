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
    plt.title("Lorenz model z versus time")

    # scatter = []
    # for n in range(len(d)):
    #     f_d = re.search(r'r-(.*?).dat', args.files[n]).group(1)
    #     plt.plot(d[n].t, d[n].z)
    # plt.show()

    phase_space_plot = [[[],[]], [[],[]]]
    for i in range(len(d[2])-1):
        if d[2].t[i]>30:
            if (d[2].x[i]*d[2].x[i+1]<0): #crossed x=0
                #interpolate
                phase_space_plot[0][0].append(0.5*(d[2].y[i]+d[2].y[i+1]))
                phase_space_plot[0][1].append(0.5*(d[2].z[i]+d[2].z[i+1]))
            if (d[2].y[i]*d[2].y[i+1]<0): #crossed y=0
                phase_space_plot[1][0].append(0.5*(d[2].z[i]+d[2].z[i+1]))
                phase_space_plot[1][1].append(0.5*(d[2].x[i]+d[2].x[i+1]))

    plt.scatter(phase_space_plot[0][0], phase_space_plot[0][1], s=0.05, c='black')
    plt.title("phase space plot: z vs y when x=0")
    plt.xlabel("y")
    plt.ylabel("z")
    plt.show()

    plt.scatter(phase_space_plot[1][0], phase_space_plot[1][1], s=0.05, c='black')
    plt.title("phase space plot: x vs z when y=0")
    plt.xlabel("z")
    plt.ylabel("x")
    plt.show()
