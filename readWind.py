# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 09:03:00 2020

@author: fritzek
"""
import pandas as pd

def readWind(sfolder,sfile,itype):
    dcswitch = {
            1: simple,
            2: swift,
            3: turbsim,
            4: mann}
        
    windFcn = dcswitch.get(itype, 'Invalid wind file type specification.')
    out = windFcn(sfolder,sfile)
    return out

def simple(sfolder,sfile):
    fwind = open(sfolder+'/'+sfile,'r')
    l1wind = fwind.readlines()
    dfwind = pd.DataFrame(columns=['time [s]','u [m/s]','v [m/s]','w [m/s]'])
    for i in range(1,len(l1wind)):
        l1line = l1wind[i].split()
        if not l1line:
            pass
        else:
            dfwind = dfwind.append(pd.Series(l1line,index=dfwind.columns),ignore_index=True)
    d2wind = dfwind.to_numpy(dtype='float')
    return d2wind

def swift(sfolder,sfile):
    return 'Swift interpreter not yet inmplemented.'

def turbsim(sfolder,sfile):
    return 'TurbSim interpreter not yet inmplemented.'

def mann(sfolder,sfile):
    return 'Mann interpreter not yet inmplemented.'