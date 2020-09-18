import pandas as pd
import matplotlib.pyplot as plt
import argparse
import re
import numpy as np

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--f', type=str,  nargs='+', help='file(s) to read')
    # parser.add_argument('--driving_freq', type=float)
    args = parser.parse_args()

    d = pd.read_csv(args.f[0])
    plt.plot(d.x, d.y)
    plt.xlabel("x (AU)")
    plt.ylabel("y (AU)")
    plt.show()
