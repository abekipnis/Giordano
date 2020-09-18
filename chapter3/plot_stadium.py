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

    a = float(re.search(r'a-(.*?)_', args.files[0]).group(1))

    d = pd.read_csv(args.files[0])
    plt.scatter(d.x, d.y, s=0.02)


    d1 = pd.read_csv(args.files[1])
    plt.scatter(d1.x, d1.y, s=0.02)

    plt.gca().set_aspect('equal')# adjustable='box')

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("trajectory for 'stadium' billiard\nalpha = %lf" %(float(a)))

    trange = np.arange(0, np.pi, 0.01)
    plt.plot(np.cos(trange), np.sin(trange)+a, color='r')#, s=0.05)
    trange = np.arange(np.pi, 2*np.pi, 0.01)
    plt.plot(np.cos(trange), np.sin(trange)-a, color='r')#, s=0.05)

    plt.plot([-1,-1],[a,-a],  color='r')
    plt.plot([1,1],[a,-a],  color='r')

    plt.show()

    # plt.plot(np.sqrt())
