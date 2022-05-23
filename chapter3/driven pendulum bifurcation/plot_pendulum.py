import pandas as pd
import matplotlib.pyplot as plt
import argparse
import re

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--files', type=str,  nargs='+', help='file(s) to read')
    parser.add_argument('--driving_freq', type=float)
    args = parser.parse_args()
    d = []

    fig, axes = plt.subplots(2, 2, sharex=True, sharey=True)
    d = [pd.read_csv(f) for f in args.files]
    for n in range(len(d)):
        phase = re.search(r'f_d-(.*?).dat', args.files[n]).group(1)
        fig.suptitle("Poincare sections at different points in the driving cycle")

        axes[n//2][n%2].scatter(d[n].theta, d[n].omega, s=0.05)
        axes[n//2][n%2].set_title("phi = %s" %(phase)) 
        axes[n//2][n%2].set_xlabel('theta')
        axes[n//2][n%2].set_ylabel('omega')
    plt.show()


    #     plt.plot(d[n].x, d[n].y)
    # plt.show()
    #
    # plt.plot(d[0].x, abs(d[0].y-d[1].y))
    # plt.yscale('log')
    # plt.show()
