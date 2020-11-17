# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 11:50:34 2020

@author: fritzek
"""
#%% TO DO
"""
include:
    structural properties
    aeroelastic functions
"""

#%% import packages
import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#from drawnow import drawnow, figure
#from tabulate import tabulate

#%% class definition
class Blade:
    def __init__(self, turbine, iBladeNum, dcinput):
        self.aeroRoot       = float(dcinput['AEROROOT'])
        self.bAxX           = np.array(dcinput['AEROPROPS']['xB'],dtype='float')
        self.bAxY           = np.array(dcinput['AEROPROPS']['yB'],dtype='float')
        self.bAxZ           = np.array(dcinput['AEROPROPS']['zB'],dtype='float')
        self.bladeRoot      = float(dcinput['BLADEROOT'])
        self.d1cl           = {}
        self.d1L            = {}
        self.d1nX           = {}
        self.d1nY           = {}
        self.d1nZ           = {}
        self.d1velIndB2Bu   = {}
        self.d1velIndB2Bv   = {}
        self.d1velIndB2Bw   = {}
        self.d1velIndW2Bu   = {}
        self.d1velIndW2Bv   = {}
        self.d1velIndW2Bw   = {}
        self.d2Fax          = {}
        self.d2Fx           = {}
        self.d2FxPerUnitL   = {}
        self.d2Fy           = {}
        self.d2FyPerUnitL   = {}
        self.d2Fz           = {}
        self.d2FzPerUnitL   = {}
        self.d2L            = {}
        self.d2velIndB2Bu   = {}
        self.d2velIndB2Bv   = {}
        self.d2velIndB2Bw   = {}
#        self.d2velIndW2Bu   = {}
#        self.d2velIndW2Bv   = {}
#        self.d2velIndW2Bw   = {}
        self.d2X            = {}
        self.d2Xcp          = {}
        self.d2Y            = {}
        self.d2Ycp          = {}
        self.d2Z            = {}
        self.d2Zcp          = {}
        self.dFaxTot        = {}
        self.dLtot          = {}
        self.Gamma          = {}
        self.GammaCor       = {}
        self.iBladeNum      = iBladeNum
        self.length         = float(dcinput['BLADELENGTH'])
        self.pitch          = float(dcinput['PITCHANGLE'])
        
    def calcNvec(self,t):
        d2AX = self.d2X[t][0:-1,1:]
        d2AY = self.d2Y[t][0:-1,1:]
        d2AZ = self.d2Z[t][0:-1,1:]
        
        d2BX = self.d2X[t][0:-1,0:-1]
        d2BY = self.d2Y[t][0:-1,0:-1]
        d2BZ = self.d2Z[t][0:-1,0:-1]
        
        d2CX = self.d2X[t][1:,0:-1]
        d2CY = self.d2Y[t][1:,0:-1]
        d2CZ = self.d2Z[t][1:,0:-1]
        
        d1DX = self.d2X[t][1:,1:]
        d1DY = self.d2Y[t][1:,1:]
        d1DZ = self.d2Z[t][1:,1:]
        
        d2nX = (d2CY-d2AY)*(d1DZ-d2BZ)-(d2CZ-d2AZ)*(d1DY-d2BY)
        d2nY = (d2CZ-d2AZ)*(d1DX-d2BX)-(d2CX-d2AX)*(d1DZ-d2BZ)
        d2nZ = (d2CX-d2AX)*(d1DY-d2BY)-(d2CY-d2AY)*(d1DX-d2BX)
        
        #%%test #TODO: delete?
        d2nNorm = np.sqrt(np.square(d2nX)+np.square(d2nY)+np.square(d2nZ))
        d2nX = d2nX/d2nNorm
        d2nY = d2nY/d2nNorm
        d2nZ = d2nZ/d2nNorm
        #%%
        
        self.d1nX[t] = np.reshape(d2nX,(1,np.size(d2nX)),order='F').T
        self.d1nY[t] = np.reshape(d2nY,(1,np.size(d2nY)),order='F').T
        self.d1nZ[t] = np.reshape(d2nZ,(1,np.size(d2nZ)),order='F').T

    def calcBladeLoads(self,sim,t,turbine,i):
        d2Gamma = self.Gamma[t].reshape([self.iElemJ,self.iElemI],order='F')
        d2GammaCor = np.hstack([d2Gamma[:,0].reshape([self.iElemJ,1]), d2Gamma[:,1:] - d2Gamma[:,0:-1]])
        
        self.GammaCor[t] = d2GammaCor.reshape([1,np.size(d2Gamma)],order='F')
        
        #spanwise distance per panel
        d2dXs = self.d2X[t][1:,0:-1]-self.d2X[t][0:-1,0:-1] 
        d2dYs = self.d2Y[t][1:,0:-1]-self.d2Y[t][0:-1,0:-1]
        d2dZs = self.d2Z[t][1:,0:-1]-self.d2Z[t][0:-1,0:-1]
        
        #chordwise distance per panel #TODO: potentially delete nect block
#        d2dXc = 0.5*(self.d2X[t][0:-1,1:]-self.d2X[t][0:-1,0:-1] + self.d2X[t][1:,1:]-self.d2X[t][1:,0:-1])
#        d2dYc = 0.5*(self.d2Y[t][0:-1,1:]-self.d2Y[t][0:-1,0:-1] + self.d2Y[t][1:,1:]-self.d2Y[t][1:,0:-1])
#        d2dZc = 0.5*(self.d2Z[t][0:-1,1:]-self.d2Z[t][0:-1,0:-1] + self.d2Z[t][1:,1:]-self.d2Z[t][1:,0:-1])
#        d2dc  = np.sqrt(np.square(d2dXc)+np.square(d2dYc)+np.square(d2dZc)) 
        
        d1dXc = 0.5*(self.d2X[t][0:-1,-1]-self.d2X[t][0:-1,0] + self.d2X[t][1:,-1]-self.d2X[t][1:,0])
        d1dYc = 0.5*(self.d2Y[t][0:-1,-1]-self.d2Y[t][0:-1,0] + self.d2Y[t][1:,-1]-self.d2Y[t][1:,0])
        d1dZc = 0.5*(self.d2Z[t][0:-1,-1]-self.d2Z[t][0:-1,0] + self.d2Z[t][1:,-1]-self.d2Z[t][1:,0])
        d1dc  = np.sqrt(np.square(d1dXc)+np.square(d1dYc)+np.square(d1dZc)) 
        
        d2eNorm = np.sqrt(np.square(d2dXs)+np.square(d2dYs)+np.square(d2dZs))

        d2e1 = d2dXs/d2eNorm
        d2e2 = d2dYs/d2eNorm
        d2e3 = d2dZs/d2eNorm
        
        d1velRot = 2*np.pi*turbine.rpm/60*np.sqrt(np.square(self.d2Ycp[t][:,int(np.floor(self.iElemI/2))])+
                                                  np.square(self.d2Zcp[t][:,int(np.floor(self.iElemI/2))]-turbine.hubHeight-turbine.zNac2hub))
        d1velRotY = np.tile(-d1velRot*np.cos(np.deg2rad(turbine.bladePos[t][i])),(self.iElemI))
        d1velRotZ = np.tile(-d1velRot*np.sin(np.deg2rad(turbine.bladePos[t][i])),(self.iElemI))
        
        d1VtotU = (sim.d2wind[0,1] + self.d1velIndB2Bu[t] + self.d1velIndW2Bu[t])
        d1VtotV = (sim.d2wind[0,2] + self.d1velIndB2Bv[t] + self.d1velIndW2Bv[t] + d1velRotY)
        d1VtotW = (sim.d2wind[0,3] + self.d1velIndB2Bw[t] + self.d1velIndW2Bw[t] + d1velRotZ)
        
        d2VtotU = np.reshape(d1VtotU,(self.iElemJ,self.iElemI),order='F')
        d2VtotV = np.reshape(d1VtotV,(self.iElemJ,self.iElemI),order='F')
        d2VtotW = np.reshape(d1VtotW,(self.iElemJ,self.iElemI),order='F')
        d2Veff  = np.sqrt(np.square(d2VtotU)+np.square(d2VtotV)+np.square(d2VtotW))

        self.d2FxPerUnitL[t] = sim.rho*(d2VtotV*d2e3 - d2VtotW*d2e2)*d2GammaCor
        self.d2FyPerUnitL[t] = sim.rho*(d2VtotW*d2e1 - d2VtotU*d2e3)*d2GammaCor
        self.d2FzPerUnitL[t] = sim.rho*(d2VtotU*d2e2 - d2VtotV*d2e1)*d2GammaCor
       
        self.d2Fx[t] = sim.rho*(d2VtotV*d2dZs - d2VtotW*d2dYs)*d2GammaCor
        self.d2Fy[t] = sim.rho*(d2VtotW*d2dXs - d2VtotU*d2dZs)*d2GammaCor
        self.d2Fz[t] = sim.rho*(d2VtotU*d2dYs - d2VtotV*d2dXs)*d2GammaCor
        
        self.d2L[t] = sim.rho*d2Veff*d2GammaCor
        self.d1L[t] = np.sum(self.d2L[t],axis=1)
        self.dLtot[t] = np.sum(np.sum(self.d2L[t]*d2eNorm,axis=1))
        
        d2cl = self.d2L[t]/(0.5*sim.rho*np.square(d2Veff))
        self.d1cl[t] = np.sum(d2cl,axis=1)/d1dc

        self.d2Fax[t] = self.d2Fx[t]
        self.dFaxTot[t] = np.sum(np.sum(self.d2Fx[t]*d2eNorm,axis=1))
        
        