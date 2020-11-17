# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 16:20:05 2020

@author: fritzek
"""
#%%packages
import os
import numpy as np
import matplotlib.pyplot as plt

os.chdir("C:/Users/fritzek/OneDrive - TNO/PhD/Code development/Python/LiftingSurfaceCodeRot")
from CSimulation import Simulation
from CTurbine import Turbine
from CPostProcessor import PostProcessor
from writeOutput import writeOutput
from readInput import readInput
from readWind import readWind

#%%input
sfolder = 'C:/Users/fritzek/OneDrive - TNO/PhD/Code development/Python/LiftingSurfaceCodeRot'
sfile = 'inputStraightNoTwist.txt'
#sairfoilFile = 'GE63.7airfoils.dat'
sairfoilFile = 'NACA4412airfoils.dat'


iElemI = 10 #chordwise discretisation
iElemJ = 40 #spanwise discretisation
#%% load data
#sfolder ='C:/Users/fritzek/OneDrive - TNO/PhD/Code development/Python'
dcinput = readInput(sfolder,sfile)

#%%Adjustments # TODO: DELETE ME
dcinput['TEND'] = '60.0'
dcinput['TIMESTEP'] = '1.0' 
dcinput['RPM'] = '0.0'
dcinput['NROFBLADES'] = '1'
dcinput['PITCHANGLE'] = '-90'
dcinput['STARTAZIMUTH'] = '-90.0'
dcinput['ELEMI'] = str(iElemI)
dcinput['ELEMJ'] = str(iElemJ)
#%% Setup simulation
sim = Simulation(dcinput,iElemI,iElemJ)
sim.d2wind = readWind(sfolder,dcinput['WINDFILENAME'],int(dcinput['WINDTYPE']))

turbine = Turbine(dcinput)
turbine.initBlades(dcinput,iElemI,iElemJ)
turbine.initWake()
turbine.discrBlade(sfolder,sairfoilFile,np.array(dcinput['AEROPROPS'],dtype='float'),iElemI,iElemJ)

#%% Run simulation
sim.run(turbine)

#%% Postprocess
tPlot = sim.tsteps-1

post = PostProcessor()
post.plotCirculation(iElemI,iElemJ,turbine,tPlot)
post.plotLift(turbine,tPlot)
post.plotCl(turbine,tPlot)
d1Mpitch = post.plotPitchMoment(turbine,tPlot)
post.plotGlobalLoads(sim,turbine,tPlot)
#post.plotVelPlane(sim,turbine,tPlot,0.,0.,0.,'X',2.)
#post.plotVelPlane(sim,turbine,tPlot,1.,0.,0.,'Y',2.)
#post.plotVelPlane(sim,turbine,tPlot,1.,0.,100.,'Z',2.)
#eng = post.plotWake(iElemI,iElemJ,turbine,tPlot)

#%%writing output
sfolderOut = sfolder+'/testOut/'
writeOutput(sfolderOut,dcinput,sim,turbine)
    
#%%delete block (only used to shortcut to non-rotating code)
d1Rz = turbine.blade[0].d1spanPos[0:-1] + 0.5*np.diff(turbine.blade[0].d1spanPos) - turbine.blade[0].bladeRoot
d1Rg = np.sum(turbine.blade[0].GammaCor[tPlot].reshape([iElemJ,iElemI],order='F'),axis=1)
d1Rc = turbine.blade[0].d1cl[tPlot]
d2out = np.vstack([d1Rz,d1Rg,d1Rc]).T
d2out1 = np.vstack([d1Rz,d1Rg,d1Rc])
np.savetxt('testStraightNoTwist.txt',d2out)
i = 0
#%%
