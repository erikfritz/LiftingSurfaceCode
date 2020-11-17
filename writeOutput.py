# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 10:16:15 2020

@author: fritzek
"""
#%% packages
import os
import numpy as np
import shutil
import json

#%% function
def writeOutput(sfolder,dcinput,sim,turbine):
    tOut = sim.tsteps - 1
    bisDir = os.path.isdir(sfolder)
    if bisDir:
        shutil.rmtree(sfolder)
        os.mkdir(sfolder)
    else:
        os.mkdir(sfolder)
        
    writeGeomFiles(sfolder,sim,turbine,tOut)
    writeWakeCirc(sfolder,turbine,tOut)
    writeBladeProps(sfolder,sim,turbine,tOut)
#    writeSimSetup(sfolder,dcinput)
    
def writeSimSetup(sfolder,dcinput):
    dctmp = dcinput
    del dctmp['AEROPROPS']
    
    ssimSetupFile = 'SimSetup.dat'
    with open(sfolder+ssimSetupFile,'w') as f:
#        print(dctmp,file=f)
        f.write(json.dumps(dctmp))
       
def writeGeomFiles(sfolder,sim,turbine,tOut):
    for i in range(0,turbine.nBlades):
        d2geomBlade = np.vstack([turbine.blade[i].d2X[tOut].reshape((sim.iElemI+1)*(sim.iElemJ+1)),
                                 turbine.blade[i].d2Y[tOut].reshape((sim.iElemI+1)*(sim.iElemJ+1)),
                                 turbine.blade[i].d2Z[tOut].reshape((sim.iElemI+1)*(sim.iElemJ+1))]).T
        sgeomBladeFile = 'B%i_bladeGeom.dat' %i
        
        d2geomWake = np.vstack([turbine.wake[i].d2X[tOut].reshape((tOut+1)*(sim.iElemJ+1)),
                                turbine.wake[i].d2Y[tOut].reshape((tOut+1)*(sim.iElemJ+1)),
                                turbine.wake[i].d2Z[tOut].reshape((tOut+1)*(sim.iElemJ+1))]).T
        sgeomWakeFile = 'B%i_wakeGeom.dat' %i
        
        np.savetxt(sfolder+sgeomBladeFile,d2geomBlade)
        np.savetxt(sfolder+sgeomWakeFile,d2geomWake)

def writeWakeCirc(sfolder,turbine,tOut):
    for i in range(0,turbine.nBlades):
        d1circWake = turbine.wake[i].Gamma[tOut]
        scircWakeFile = 'B%i_wakeCirc.dat' %i
        
        np.savetxt(sfolder+scircWakeFile,d1circWake)
        
def writeBladeProps(sfolder,sim,turbine,tOut):
    for i in range(0,turbine.nBlades):
        d2bladeProps = np.vstack([turbine.blade[i].d2Xcp[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),            #1
                                  turbine.blade[i].d2Ycp[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),            #2
                                  turbine.blade[i].d2Zcp[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),            #3
                                  turbine.blade[i].d1nX[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),             #4
                                  turbine.blade[i].d1nY[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),             #5
                                  turbine.blade[i].d1nZ[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),             #6
                                  turbine.blade[i].Gamma[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),            #7
                                  turbine.blade[i].GammaCor[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),         #8
                                  turbine.blade[i].d2Fx[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),             #9
                                  turbine.blade[i].d2Fy[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),             #10
                                  turbine.blade[i].d2Fz[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),             #11
                                  turbine.blade[i].d2FxPerUnitL[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),     #12
                                  turbine.blade[i].d2FyPerUnitL[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),     #13
                                  turbine.blade[i].d2FzPerUnitL[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),     #14
                                  turbine.blade[i].d2L[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),              #15
                                  turbine.blade[i].d1velIndB2Bu[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),     #16
                                  turbine.blade[i].d1velIndB2Bv[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),     #17
                                  turbine.blade[i].d1velIndB2Bw[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),     #18
                                  turbine.blade[i].d1velIndW2Bu[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),     #19
                                  turbine.blade[i].d1velIndW2Bv[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F'),     #20
                                  turbine.blade[i].d1velIndW2Bw[tOut].reshape((sim.iElemI)*(sim.iElemJ),order='F')]).T  #21
        sbladePropsFile = 'B%i_bladeProps.dat' %i
        
        np.savetxt(sfolder+sbladePropsFile,d2bladeProps)