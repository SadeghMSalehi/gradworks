__author__ = 'joohwi'

import numpy as np
import matplotlib.pyplot as plt


def two_scale(t, s1, s2, fname = ""):
    fig, ax1 = plt.subplots()

    ax1.plot(t, s1, 'b-')
    ax1.set_xlabel('time (s)')
    # Make the y-axis label and tick labels match the line color.
    ax1.set_ylabel('exp', color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')

    ax2 = ax1.twinx()
    ax2.plot(t, s2, 'r.')
    ax2.set_ylabel('sin', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')

    if fname == "":
        plt.show()
    else:
        plt.savefig(fname)



t = np.arange(0.01, 10.0, 0.01)
s1 = np.exp(t)
s2 = np.sin(2*np.pi*t)
two_scale(t, s1, s2, "/Prime/Thesis-Data/SimilarityMetric/twoscale_example.pdf")