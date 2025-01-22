"""
This python file includes the packages needed to be imported for the data loading , computations,
and analyses
"""

import os
from itertools import combinations
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.io import loadmat
from abc import ABC, abstractmethod
from scipy.stats import multivariate_normal
from matplotlib.ticker import MaxNLocator



from scipy.ndimage import gaussian_filter1d
from matplotlib.patches import Patch
from scipy.stats import ttest_ind, shapiro