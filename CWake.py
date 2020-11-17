# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 13:48:15 2020

@author: fritzek
"""

#%% TO DO
"""
include:
"""

#%% import packages
import numpy as np
#import time
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

#%% class definition
class Wake:
    def __init__(self, turbine, iWake):
        self.iWake = iWake
        self.d1velIndB2Wu = {}
        self.d1velIndB2Wv = {}
        self.d1velIndB2Ww = {}
        self.d1velIndW2Wu = {}
        self.d1velIndW2Wv = {}
        self.d1velIndW2Ww = {}
        self.d2X   = {} 
        self.d2Y   = {} 
        self.d2Z   = {} 
        self.Gamma = {}
