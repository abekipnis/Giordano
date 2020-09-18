import pandas as pd
import matplotlib.pyplot as plt
import argparse
import re
import numpy as np
import scipy.optimize

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--files', type=str,  nargs='+', help='file(s) to read')
    # parser.add_argument('--driving_freq', type=float)
    args = parser.parse_args()


    d = pd.read_csv(args.files[0], names=['t','ds'], header=None)
    plt.scatter(d.t, d.ds, s=0.02)


    fit, fit2 = scipy.optimize.curve_fit(lambda t,a,b,c: a+np.exp(b*t),  d.t,  d.ds)

    f = np.polyfit(d.t[1:-4], np.log(d.ds[1:-4]), 1)
    plt.plot(d.t, np.exp(f[1])*np.exp(f[0]*d.t), c='red', label=r'fit, $\lambda$=%lf' %(f[0]))


    plt.xlabel("time (s)")
    plt.ylabel("average difference between trajectories s(t)")
    plt.yscale("log")
    plt.legend()
    alpha = float(re.search(r'a-(.*?).dat', args.files[0]).group(1))

    plt.savefig("Stadium Lyapunov exponent alpha=%.10lf.png" %(alpha))
    # plt.show()
    print("%lf" %(f[0]))
