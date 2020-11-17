# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 11:44:07 2020

@author: fritzek
"""

#%% TO DO
"""
include:
"""

#%% import packages
import numpy as np
import time
from vortexRing import vortexRingVec

#from vortexRing import vortexRingVec
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

#%% class definition
class Simulation:
    def __init__(self, dcinput,iElemI,iElemJ):
        self.t0             = float(dcinput['TBEGIN'])
        self.tend           = float(dcinput['TEND'])
        self.delT           = float(dcinput['TIMESTEP'])
        self.t              = np.arange(self.t0,self.tend,self.delT)
        self.tsteps         = len(self.t)
        self.rho            = float(dcinput['AIRDENSITY'])
        self.wallT          = np.zeros_like(self.t)
        self.wakePoints     = int(dcinput['STREAMWISEWAKEPOINTS'])
        self.freeWakePoints = int(dcinput['FREESTRMWISEWAKEPOINTS'])
        self.iElemI         = iElemI
        self.iElemJ         = iElemJ

    def run(self,turbine):
        tic = time.time()
        print('Progress   Simulated time   Wall time   No. of vortex elements')
        
        self.tstep0(turbine)
        
        for t in range(1,self.tsteps):
            self.tstep(t,turbine)
            self.wallT[t] = time.time() - tic
            
            print(str("{:.2f}".format(round(100*(t+1)/self.tsteps,2)))+'%\t\t'+
                  str("{:.2f}".format(round(self.t[t],3)))+'s\t\t'+
                  time.strftime('%H:%M:%S',time.gmtime(self.wallT[t]))+'\t'+
                  str(turbine.nBlades*(self.iElemI*self.iElemJ+(t-1)*self.iElemJ)))

    def tstep0(self,turbine):
        for i in range(0,turbine.nBlades):
            turbine.wake[i].d2X[0] = turbine.blade[i].d2X[0][:,-1]
            turbine.wake[i].d2Y[0] = turbine.blade[i].d2Y[0][:,-1]
            turbine.wake[i].d2Z[0] = turbine.blade[i].d2Z[0][:,-1]
            
            turbine.blade[i].calcNvec(0)
            
        for i in range(0,turbine.nBlades):
            for j in range(0,turbine.nBlades):
                tmpU,tmpV,tmpW,_,_,_ = vortexRingVec(turbine.blade[i].d2Xcp[0],turbine.blade[i].d2Ycp[0],turbine.blade[i].d2Zcp[0],
                                                     turbine.blade[j].d2X[0],turbine.blade[j].d2Y[0],turbine.blade[j].d2Z[0],
                                                     np.ones(np.size(turbine.blade[i].d2Xcp[0])))
                if j == 0:
                    turbine.blade[i].d2velIndB2Bu[0] = tmpU
                    turbine.blade[i].d2velIndB2Bv[0] = tmpV
                    turbine.blade[i].d2velIndB2Bw[0] = tmpW
                else:
                    turbine.blade[i].d2velIndB2Bu[0] = turbine.blade[i].d2velIndB2Bu[0] + tmpU
                    turbine.blade[i].d2velIndB2Bv[0] = turbine.blade[i].d2velIndB2Bv[0] + tmpV
                    turbine.blade[i].d2velIndB2Bw[0] = turbine.blade[i].d2velIndB2Bw[0] + tmpW
                    
            d2aIJ = (np.tile(turbine.blade[i].d1nX[0],(1,self.iElemI*self.iElemJ))*turbine.blade[i].d2velIndB2Bu[0] + 
                        np.tile(turbine.blade[i].d1nY[0],(1,self.iElemI*self.iElemJ))*turbine.blade[i].d2velIndB2Bv[0] + 
                        np.tile(turbine.blade[i].d1nZ[0],(1,self.iElemI*self.iElemJ))*turbine.blade[i].d2velIndB2Bw[0])

            d1velRot = 2*np.pi*turbine.rpm/60*np.sqrt(np.square(turbine.blade[i].d2Ycp[0][:,int(np.floor(self.iElemI/2))])+
                                                      np.square(turbine.blade[i].d2Zcp[0][:,int(np.floor(self.iElemI/2))]-turbine.hubHeight-turbine.zNac2hub))
            d1velRotY = np.tile(d1velRot*np.cos(np.deg2rad(turbine.bladePos[0][i])),(self.iElemI))
            d1velRotZ = np.tile(d1velRot*np.sin(np.deg2rad(turbine.bladePos[0][i])),(self.iElemI))

            d1RHS = -(self.d2wind[0,1]*turbine.blade[i].d1nX[0]+ 
                      self.d2wind[0,2]*turbine.blade[i].d1nY[0]+ 
                      self.d2wind[0,3]*turbine.blade[i].d1nZ[0]+
                      np.reshape(d1velRotY,(np.size(d1velRotY),1))*turbine.blade[i].d1nY[0]+
                      np.reshape(d1velRotZ,(np.size(d1velRotZ),1))*turbine.blade[i].d1nZ[0]) #TODO: needs adjustment
            
            turbine.blade[i].Gamma[0] = np.linalg.solve(d2aIJ, d1RHS)

            d2Gamma = turbine.blade[i].Gamma[0].reshape([self.iElemJ,self.iElemI],order='F')
            d2GammaCor = np.hstack([d2Gamma[:,0].reshape([self.iElemJ,1]), d2Gamma[:,1:] - d2Gamma[:,0:-1]])
            turbine.blade[i].GammaCor[0] = d2GammaCor.reshape([1,np.size(d2Gamma)],order='F')
            
    def tstep(self,t,turbine):
        turbine.rotate(t,self.delT)
        #W2W
        if t == 1:
            for i in range(0,turbine.nBlades):
                turbine.wake[i].d1velIndW2Wu[t] = np.zeros_like(turbine.blade[i].d2X[t][:,-1])
                turbine.wake[i].d1velIndW2Wv[t] = np.zeros_like(turbine.blade[i].d2Y[t][:,-1])
                turbine.wake[i].d1velIndW2Ww[t] = np.zeros_like(turbine.blade[i].d2Z[t][:,-1])
                
                turbine.wake[i].Gamma[t] = turbine.blade[i].Gamma[t-1][-self.iElemJ:]
        else:
            for i in range(0,turbine.nBlades):
                for j in range(0,turbine.nBlades):
                    tmpU,tmpV,tmpW,_,_,_ = vortexRingVec(turbine.wake[i].d2X[t-1],turbine.wake[i].d2Y[t-1],turbine.wake[i].d2Z[t-1],
                                                         turbine.wake[j].d2X[t-1],turbine.wake[j].d2Y[t-1],turbine.wake[j].d2Z[t-1],
                                                         turbine.wake[j].Gamma[t-1])
                    if j == 0:
                        d2velIndW2Wu = tmpU
                        d2velIndW2Wv = tmpV
                        d2velIndW2Ww = tmpW
                    else:
                        d2velIndW2Wu = d2velIndW2Wu + tmpU
                        d2velIndW2Wv = d2velIndW2Wv + tmpV
                        d2velIndW2Ww = d2velIndW2Ww + tmpW

                turbine.wake[i].d1velIndW2Wu[t] = np.sum(d2velIndW2Wu,axis=1)
                turbine.wake[i].d1velIndW2Wv[t] = np.sum(d2velIndW2Wv,axis=1)
                turbine.wake[i].d1velIndW2Ww[t] = np.sum(d2velIndW2Ww,axis=1)

                turbine.wake[i].Gamma[t] = np.vstack([turbine.blade[i].Gamma[t-1][-self.iElemJ:],turbine.wake[i].Gamma[t-1]])
                
        #B2W
        for i in range(0,turbine.nBlades):
            for j in range(0,turbine.nBlades):
                tmpU,tmpV,tmpW,_,_,_ = vortexRingVec(turbine.wake[i].d2X[t-1],turbine.wake[i].d2Y[t-1],turbine.wake[i].d2Z[t-1],
                                                     turbine.blade[j].d2X[t-1],turbine.blade[j].d2Y[t-1],turbine.blade[j].d2Z[t-1],
                                                     turbine.blade[j].Gamma[t-1])
 
                if j == 0:
                    d2velIndB2Wu = tmpU
                    d2velIndB2Wv = tmpV
                    d2velIndB2Ww = tmpW
                else:
                    d2velIndB2Wu = d2velIndB2Wu + tmpU
                    d2velIndB2Wv = d2velIndB2Wv + tmpV
                    d2velIndB2Ww = d2velIndB2Ww + tmpW
                    
            turbine.wake[i].d1velIndB2Wu[t] = np.sum(d2velIndB2Wu,axis=1)
            turbine.wake[i].d1velIndB2Wv[t] = np.sum(d2velIndB2Wv,axis=1)
            turbine.wake[i].d1velIndB2Ww[t] = np.sum(d2velIndB2Ww,axis=1)

        #Convect wake
        for i in range(0,turbine.nBlades):
            d2wakeConvX = np.reshape((self.d2wind[0,1]*np.ones_like(turbine.wake[i].d1velIndB2Wu[t])+
                                      turbine.wake[i].d1velIndW2Wu[t]+
                                      turbine.wake[i].d1velIndB2Wu[t])*self.delT,[self.iElemJ+1,t],order='F') #TODO: needs adjustment
            
            d2wakeConvY = np.reshape((self.d2wind[0,2]*np.ones_like(turbine.wake[i].d1velIndB2Wv[t])+
                                      turbine.wake[i].d1velIndW2Wv[t]+
                                      turbine.wake[i].d1velIndB2Wv[t])*self.delT,[self.iElemJ+1,t],order='F')
            
            d2wakeConvZ = np.reshape((self.d2wind[0,3]*np.ones_like(turbine.wake[i].d1velIndB2Ww[t])+
                                      turbine.wake[i].d1velIndW2Ww[t]+
                                      turbine.wake[i].d1velIndB2Ww[t])*self.delT,[self.iElemJ+1,t],order='F')
         
            turbine.wake[i].d2X[t] = np.append([turbine.blade[i].d2X[t][:,-1]],turbine.wake[i].d2X[t-1].T+d2wakeConvX.T,axis=0).T
            turbine.wake[i].d2Y[t] = np.append([turbine.blade[i].d2Y[t][:,-1]],turbine.wake[i].d2Y[t-1].T+d2wakeConvY.T,axis=0).T
            turbine.wake[i].d2Z[t] = np.append([turbine.blade[i].d2Z[t][:,-1]],turbine.wake[i].d2Z[t-1].T+d2wakeConvZ.T,axis=0).T
            
        #W2B
        for i in range(0,turbine.nBlades):
            for j in range(0,turbine.nBlades):
                tmpU,tmpV,tmpW,_,_,_ = vortexRingVec(turbine.blade[i].d2Xcp[t],turbine.blade[i].d2Ycp[t],turbine.blade[i].d2Zcp[t],
                                                     turbine.wake[j].d2X[t],turbine.wake[j].d2Y[t],turbine.wake[j].d2Z[t],
                                                     turbine.wake[j].Gamma[t])
 
                if j == 0:
                    d2velIndW2Bu = tmpU
                    d2velIndW2Bv = tmpV
                    d2velIndW2Bw = tmpW
                else:
                    d2velIndW2Bu = d2velIndW2Bu + tmpU
                    d2velIndW2Bv = d2velIndW2Bv + tmpV
                    d2velIndW2Bw = d2velIndW2Bw + tmpW
                    
            turbine.blade[i].d1velIndW2Bu[t] = np.sum(d2velIndW2Bu,axis=1)
            turbine.blade[i].d1velIndW2Bv[t] = np.sum(d2velIndW2Bv,axis=1)
            turbine.blade[i].d1velIndW2Bw[t] = np.sum(d2velIndW2Bw,axis=1)
        
        #B2B
        for i in range(0,turbine.nBlades):
            for j in range(0,turbine.nBlades):
                tmpU,tmpV,tmpW,_,_,_ = vortexRingVec(turbine.blade[i].d2Xcp[t],turbine.blade[i].d2Ycp[t],turbine.blade[i].d2Zcp[t],
                                                     turbine.blade[j].d2X[t],turbine.blade[j].d2Y[t],turbine.blade[j].d2Z[t],
                                                     np.ones(np.size(turbine.blade[i].d2Xcp[t])))
                #%% test
#                d2testX = np.vstack((turbine.blade[i].d2X[t].T,turbine.wake[i].d2X[t][:,1])).T
#                d2testY = np.vstack((turbine.blade[i].d2Y[t].T,turbine.wake[i].d2Y[t][:,1])).T
#                d2testZ = np.vstack((turbine.blade[i].d2Z[t].T,turbine.wake[i].d2Z[t][:,1])).T
#                tmpU,tmpV,tmpW,_,_,_ = vortexRingVec(turbine.blade[i].d2Xcp[t],turbine.blade[i].d2Ycp[t],turbine.blade[i].d2Zcp[t],
#                                                     d2testX,d2testY,d2testZ,
#                                                     np.ones((self.iElemJ*(self.iElemI+1))))
#                tmpUwake = np.sum(tmpU[:,-self.iElemJ:],axis=1)
#                tmpVwake = np.sum(tmpV[:,-self.iElemJ:],axis=1)
#                tmpWwake = np.sum(tmpW[:,-self.iElemJ:],axis=1)
#                tmpU = tmpU[:,0:-self.iElemJ]
#                tmpV = tmpV[:,0:-self.iElemJ]
#                tmpW = tmpW[:,0:-self.iElemJ]
#                tmpU[:,-1] = tmpU[:,-1] + tmpUwake
#                tmpV[:,-1] = tmpV[:,-1] + tmpVwake
#                tmpW[:,-1] = tmpW[:,-1] + tmpWwake
                #%%
 
                if j == 0:
                    turbine.blade[i].d2velIndB2Bu[t] = tmpU
                    turbine.blade[i].d2velIndB2Bv[t] = tmpV
                    turbine.blade[i].d2velIndB2Bw[t] = tmpW
                else:
                    turbine.blade[i].d2velIndB2Bu[t] = turbine.blade[i].d2velIndB2Bu[t] + tmpU
                    turbine.blade[i].d2velIndB2Bv[t] = turbine.blade[i].d2velIndB2Bv[t] + tmpV
                    turbine.blade[i].d2velIndB2Bw[t] = turbine.blade[i].d2velIndB2Bw[t] + tmpW
                    
            turbine.blade[i].d1velIndB2Bu[t] = np.sum(turbine.blade[i].d2velIndB2Bu[t],axis=1)
            turbine.blade[i].d1velIndB2Bv[t] = np.sum(turbine.blade[i].d2velIndB2Bv[t],axis=1)
            turbine.blade[i].d1velIndB2Bw[t] = np.sum(turbine.blade[i].d2velIndB2Bw[t],axis=1)
        
        #Circulation
        for i in range(0,turbine.nBlades):
            turbine.blade[i].calcNvec(t)
            
            d2aIJ = (np.tile(turbine.blade[i].d1nX[t],(1,self.iElemI*self.iElemJ))*turbine.blade[i].d2velIndB2Bu[t] + 
                     np.tile(turbine.blade[i].d1nY[t],(1,self.iElemI*self.iElemJ))*turbine.blade[i].d2velIndB2Bv[t] + 
                     np.tile(turbine.blade[i].d1nZ[t],(1,self.iElemI*self.iElemJ))*turbine.blade[i].d2velIndB2Bw[t])
            
            #TODO: wind velocity needs adjustment
            d1velRot = 2*np.pi*turbine.rpm/60*np.sqrt(np.square(turbine.blade[i].d2Ycp[t][:,int(np.floor(self.iElemI/2))])+
                                                      np.square(turbine.blade[i].d2Zcp[t][:,int(np.floor(self.iElemI/2))]-turbine.hubHeight-turbine.zNac2hub))
            d1velRotY = np.tile(d1velRot*np.cos(np.deg2rad(turbine.bladePos[t][i])),(self.iElemI))
            d1velRotZ = np.tile(d1velRot*np.sin(np.deg2rad(turbine.bladePos[t][i])),(self.iElemI))
            #TODO: check whether reshape is even necessary (might already have the correct shape)
            d1RHS = -(self.d2wind[0,1]*turbine.blade[i].d1nX[t]+  
                      self.d2wind[0,2]*turbine.blade[i].d1nY[t]+ 
                      self.d2wind[0,3]*turbine.blade[i].d1nZ[t]+
                      np.reshape(d1velRotY,(np.size(d1velRotY),1))*turbine.blade[i].d1nY[t]+
                      np.reshape(d1velRotZ,(np.size(d1velRotZ),1))*turbine.blade[i].d1nZ[t]+
                      np.reshape(turbine.blade[i].d1velIndW2Bu[t],(np.size(turbine.blade[i].d1velIndW2Bu[t]),1))*turbine.blade[i].d1nX[t]+
                      np.reshape(turbine.blade[i].d1velIndW2Bv[t],(np.size(turbine.blade[i].d1velIndW2Bv[t]),1))*turbine.blade[i].d1nY[t]+
                      np.reshape(turbine.blade[i].d1velIndW2Bw[t],(np.size(turbine.blade[i].d1velIndW2Bw[t]),1))*turbine.blade[i].d1nZ[t])

            turbine.blade[i].Gamma[t] = np.linalg.solve(d2aIJ, d1RHS)
            
        self.calcLoads(t,turbine)
            
    def calcLoads(self,t,turbine):
        for i in range(0,turbine.nBlades):
            turbine.blade[i].calcBladeLoads(self,t,turbine,i)
            
        turbine.calcTurbineLoads(t)
            