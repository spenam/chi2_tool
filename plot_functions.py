import boost_histogram as bh
import matplotlib.pyplot as plt
import numpy as np

def plothist(h, ax, **kwargs):
    return ax.bar(*h.axes.centers, h.values(), width=h.axes.widths[0], alpha=0.6, **kwargs)


def plothist2d(h, ax, **kwargs):
    return ax.pcolormesh(*h.axes.edges.T, h.values().T, **kwargs)