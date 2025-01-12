from src._imports import *


def set_style():
    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({'figure.titlesize': 26, 'figure.titleweight': 'bold', 'axes.titlesize': 22,
                         'axes.titleweight': "bold", 'axes.labelsize': 24,'axes.labelweight': 'bold',
                         "ytick.labelsize": 28, "xtick.labelsize": 28, 'legend.fontsize': 24,'font.family': 'Candara,  Arial'})
    plt.rcParams.update({'figure.figsize': (8, 6)})
    plt.rcParams.update({'savefig.dpi': 300})


