# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 10:25:21 2020

@author: fritzek
"""
#%% packages
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

import os
os.chdir("C:/Users/fritzek/OneDrive - TNO/PhD/Code development/Python/LiftingSurfaceCodeRot")
from vortexRing import vortexRingVec, vortexLineVec
#from readInput import readInput
#from readWind import readWind
import geomParser
os.chdir("C:/Users/fritzek/OneDrive - TNO/PhD/Code development/Python/LiftingSurfaceCodeRot")
#import time

#%% input
dalphaStart = -5.
dalphaEnd = 20.
dalphaStep = 0.5
dsweep = 0.

dconvCrit = 1e-5

dRe = 1000000
drho = 1.225
dmu = 1.802e-5

# geometry
saerofoil = 'C:/Users/fritzek/OneDrive - TNO/PhD/Code development/Python/Airfoils/NACA0012.dat'
dspan = 50.
dchord = 1.
d1wingStart = np.array([0.,0.,-dspan/2])
d1wingEnd = np.array([0.,0.,dspan/2])
iElemI = 5
iElemJ = 10

dQ = dRe*dmu/(drho*dchord)
dTimeStep = 10*dchord/dQ

#%% geometrical calculations
# wing
d1zLE = d1zTE = np.linspace(d1wingStart[2],d1wingEnd[2],100)
d1chord = dchord*np.ones_like(d1zLE) #straight chord
#d1chord = np.sqrt(1-np.square(d1zLE)/(np.square(dspan/2))) #elliptical chord
d1xLE = -0.25*d1chord
d1yLE = np.zeros_like(d1xLE)
d1xTE = 0.75*d1chord
d1yTE = np.zeros_like(d1xTE)

# airfoil
d1data = np.fromfile(saerofoil,dtype=float,sep=' ')
d1airfoilX = d1data[::2]
d1airfoilY = d1data[1::2]
iLEidx = np.argmin(d1airfoilX)
dLEval = d1airfoilX[iLEidx]
d1dataSSx = np.flip(d1airfoilX[0:iLEidx+1],0)
d1dataSSy = np.flip(d1airfoilY[0:iLEidx+1],0)
d1dataPSx = d1airfoilX[iLEidx:-1]
d1dataPSy = d1airfoilY[iLEidx:-1]

d2airfoilSS = np.vstack((np.linspace(0,1,101),np.interp(np.linspace(0,1,101),d1dataSSx,d1dataSSy)))
d2airfoilPS = np.vstack((np.linspace(0,1,101),np.interp(np.linspace(0,1,101),d1dataPSx,d1dataPSy)))
d2cambLine = np.vstack((np.linspace(0,1,101),d2airfoilPS[1,:]+0.5*(d2airfoilSS[1,:]-d2airfoilPS[1,:])))
d1cambLineGrid = np.interp(np.linspace(0,1,iElemI+1),d2cambLine[0,:],d2cambLine[1,:])

# panelling
d1spanPos = np.cos(np.linspace(np.pi,0,iElemJ+1))*dspan/2
d2gridX = np.zeros([iElemJ+1,iElemI+1])
d2gridY = np.zeros([iElemJ+1,iElemI+1])
d2gridZ = np.zeros([iElemJ+1,iElemI+1])
d2locLE = np.zeros([iElemJ+1,3])
d2locTE = np.zeros([iElemJ+1,3])
d1locChord = np.zeros(iElemJ+1)

for j in range(0,iElemJ+1):
    d2locLE[j,:] = np.array((np.interp(d1spanPos[j],d1zLE,d1xLE), np.interp(d1spanPos[j],d1zLE,d1yLE), d1spanPos[j]))
    d2locTE[j,:] = np.array((np.interp(d1spanPos[j],d1zTE,d1xTE), np.interp(d1spanPos[j],d1zTE,d1yTE), d1spanPos[j]))
    d1locChord[j] = np.sqrt(np.sum(np.square(d2locLE[j,:]-d2locTE[j,:])))
    d2locChordLine = np.linspace(d2locLE[j,:],d2locTE[j,:],iElemI+1)
    
    d2gridX[j,:] = d2locChordLine[:,0]+d1cambLineGrid*d1locChord[j]
    d2gridY[j,:] = d2locChordLine[:,1]+d1cambLineGrid*d1locChord[j]
    d2gridZ[j,:] = np.ones(iElemI+1)*d1spanPos[j]
    
d2gridY = np.nan_to_num(d2gridY)

#%%create meshgrid for VRs and CPs
#vortex rings
d2VRgridX,d2VRgridY,d2VRgridZ = geomParser.vrGrid(d2gridX,d2gridY,d2gridZ)
d2nX,d2nY,d2nZ = geomParser.nVec(d2VRgridX,d2VRgridY,d2VRgridZ)
d1nX = np.reshape(d2nX,(1,np.size(d2nX)),order='F').T
d1nY = np.reshape(d2nY,(1,np.size(d2nY)),order='F').T
d1nZ = np.reshape(d2nZ,(1,np.size(d2nZ)),order='F').T
#collocation points
d2CPgridX,d2CPgridY,d2CPgridZ=geomParser.cpGrid(d2gridX,d2gridY,d2gridZ)

#%%simulation
d1alpha = np.arange(dalphaStart,dalphaEnd+dalphaStep,dalphaStep)

d1Gamma = np.zeros_like(d1alpha)
d1cl = np.zeros_like(d1alpha)
d1cd = np.zeros_like(d1alpha)

for i in range(0,len(d1alpha)):
    dres1 = 1
    dres2 = 1
    dres3 = 1
    
    d1Q = np.array([np.cos(d1alpha[i]*np.pi/180)*np.cos(dsweep*np.pi/180)*dQ,
                    np.sin(d1alpha[i]*np.pi/180)*dQ,
                    np.cos(d1alpha[i]*np.pi/180)*np.sin(dsweep*np.pi/180)*dQ])
    
    iiter = 0
    
    st1bound = {}
    st1wake = {}
    t = 0
    
    print(' AoA\tIter\t Cl\t\t Cd\t\tRes')

    while (dres1 > dconvCrit) or (dres2 > dconvCrit) or (dres3 > dconvCrit):
        st1bound[t] = {}
        st1wake[t] = {}
    
        #influence coefficients bound vortices on themselves
        d2VelIndUB2B,d2VelIndVB2B,d2VelIndWB2B,_,_,_ = vortexRingVec(d2CPgridX,d2CPgridY,d2CPgridZ,
                                                                     d2VRgridX,d2VRgridY,d2VRgridZ,
                                                                     np.ones(np.size(d2CPgridX)))
        d1TEX = d2VRgridX[:,-1]
        d1TEY = d2VRgridY[:,-1]
        d1TEZ = d2VRgridZ[:,-1]
        #%%
#TODO: rodo FW2B (First wake row to bound circulation, fulfillment of Kutta condition at TE)
#        if t > 0:
#            d1CPgridX = d2CPgridX.reshape((1,np.size(d2CPgridX)),order='F')
#            d1CPgridY = d2CPgridY.reshape((1,np.size(d2CPgridY)),order='F')
#            d1CPgridZ = d2CPgridZ.reshape((1,np.size(d2CPgridZ)),order='F')
#            
#            d1TE1X = d1TEX[0:-1].reshape((1,iElemJ))
#            d1TE1Y = d1TEY[0:-1].reshape((1,iElemJ))
#            d1TE1Z = d1TEZ[0:-1].reshape((1,iElemJ))
#
#            d1TE2X = d1TEX[1:].reshape((1,iElemJ))
#            d1TE2Y = d1TEY[1:].reshape((1,iElemJ))
#            d1TE2Z = d1TEZ[1:].reshape((1,iElemJ))
#            
#            d2VelIndUFW2B,d2VelIndVFW2B,d2VelIndWFW2B = vortexLineVec(d1CPgridX,d1CPgridY,d1CPgridZ,
#                                                                             d1TE2X,d1TE2Y,d1TE2Z,
#                                                                             d1TE1X,d1TE1Y,d1TE1Z,
#                                                                             np.ones(iElemJ))
#            
#            d2VelIndUB2B[:,-iElemJ:] = d2VelIndUB2B[:,-iElemJ:] + d2VelIndUFW2B
#            d2VelIndVB2B[:,-iElemJ:] = d2VelIndVB2B[:,-iElemJ:] + d2VelIndVFW2B
#            d2VelIndWB2B[:,-iElemJ:] = d2VelIndWB2B[:,-iElemJ:] + d2VelIndWFW2B
#%%
            
        d2aIJ = np.tile(d1nX,(1,iElemI*iElemJ))*d2VelIndUB2B + np.tile(d1nY,(1,iElemI*iElemJ))*d2VelIndVB2B + np.tile(d1nZ,(1,iElemI*iElemJ))*d2VelIndWB2B
                    
        if t==0:
            st1wake[t]['X'] = d2VRgridX[:,-1]
            st1wake[t]['Y'] = d2VRgridY[:,-1]
            st1wake[t]['Z'] = d2VRgridZ[:,-1]
            st1wake[t]['Gamma'] = []
            
            d1RHS = -(d1Q[0]*d1nX + d1Q[1]*d1nY + d1Q[2]*d1nZ)
            st1bound[t]['Gamma'] = np.linalg.solve(d2aIJ, d1RHS)
    
        else: 
            #W2W
            if np.size(st1wake[t-1]['X']) > iElemJ+1: #wake consists of multiple rows of vortex rings
                d2velIndUW2W,d2velIndVW2W,d2velIndWW2W,_,_,_= vortexRingVec(st1wake[t-1]['X'],st1wake[t-1]['Y'],st1wake[t-1]['Z'],
                                                                            st1wake[t-1]['X'],st1wake[t-1]['Y'],st1wake[t-1]['Z'],
                                                                            st1wake[t-1]['Gamma'])
                d1velIndUW2W = np.sum(d2velIndUW2W,axis=1)
                d1velIndVW2W = np.sum(d2velIndVW2W,axis=1)
                d1velIndWW2W = np.sum(d2velIndWW2W,axis=1)
                
                st1wake[t]['Gamma'] = np.vstack([st1bound[t-1]['Gamma'][-iElemJ:],st1wake[t-1]['Gamma']])
                
            else: #wake consists of only one row of vortex rings
                d1velIndUW2W = np.zeros_like(d1TEX)
                d1velIndVW2W = np.zeros_like(d1TEY)
                d1velIndWW2W = np.zeros_like(d1TEZ)
                
                st1wake[t]['Gamma'] = st1bound[t-1]['Gamma'][-iElemJ:]
    
            #B2W
            d2velIndUB2W,d2velIndVB2W,d2velIndWB2W,_,_,_= vortexRingVec(st1wake[t-1]['X'],st1wake[t-1]['Y'],st1wake[t-1]['Z'],
                                                                        d2VRgridX,d2VRgridY,d2VRgridZ,
                                                                        st1bound[t-1]['Gamma'])
            d1velIndUB2W = np.sum(d2velIndUB2W,axis=1)
            d1velIndVB2W = np.sum(d2velIndVB2W,axis=1)
            d1velIndWB2W = np.sum(d2velIndWB2W,axis=1)
            
            #Convect wake
            d2wakeConvX = np.reshape((d1Q[0]*np.ones_like(d1velIndUB2W)+d1velIndUW2W+d1velIndUB2W)*dTimeStep,[iElemJ+1,t],order='F')
            d2wakeConvY = np.reshape((d1Q[1]*np.ones_like(d1velIndVB2W)+d1velIndVW2W+d1velIndVB2W)*dTimeStep,[iElemJ+1,t],order='F')
            d2wakeConvZ = np.reshape((d1Q[2]*np.ones_like(d1velIndWB2W)+d1velIndWW2W+d1velIndWB2W)*dTimeStep,[iElemJ+1,t],order='F')
           
            st1wake[t]['X'] = np.append([d1TEX],st1wake[t-1]['X'].T+d2wakeConvX.T,axis=0).T
            st1wake[t]['Y'] = np.append([d1TEY],st1wake[t-1]['Y'].T+d2wakeConvY.T,axis=0).T
            st1wake[t]['Z'] = np.append([d1TEZ],st1wake[t-1]['Z'].T+d2wakeConvZ.T,axis=0).T
            
            #W2B
            d2velIndUW2B,d2velIndVW2B,d2velIndWW2B,d2velDWU,d2velDWV,d2velDWW = vortexRingVec(d2CPgridX,d2CPgridY,d2CPgridZ,
                                                                                              st1wake[t]['X'],st1wake[t]['Y'],st1wake[t]['Z'],
                                                                                              st1wake[t]['Gamma'])
            #%% TODO: test whether first row of wake has to be adjusted in the same way as the TE of the bound vortices
#            d2velIndUTE2W,d2velIndVTE2W,d2velIndWTE2W = vortexLineVec(d1CPgridX,d1CPgridY,d1CPgridZ,
#                                                                      d1TE1X,d1TE1Y,d1TE1Z,
#                                                                      d1TE2X,d1TE2Y,d1TE2Z,
#                                                                      st1bound[t-1]['Gamma'][-iElemJ:])
#            d2velIndUW2B[:,0:iElemJ] = d2velIndUW2B[:,0:iElemJ] + d2velIndUTE2W
#            d2velIndVW2B[:,0:iElemJ] = d2velIndVW2B[:,0:iElemJ] + d2velIndVTE2W
#            d2velIndWW2B[:,0:iElemJ] = d2velIndWW2B[:,0:iElemJ] + d2velIndWTE2W
            #%%
            
            d1velIndUW2B = np.sum(d2velIndUW2B,axis=1)
            d1velIndVW2B = np.sum(d2velIndVW2B,axis=1)
            d1velIndWW2B = np.sum(d2velIndWW2B,axis=1)
            
            #circulation
            d1RHS = -(d1Q[0]*d1nX + d1Q[1]*d1nY + d1Q[2]*d1nZ + 
                     np.reshape(d1velIndUW2B,(np.size(d1velIndUW2B),1))*d1nX + 
                     np.reshape(d1velIndVW2B,(np.size(d1velIndVW2B),1))*d1nY + 
                     np.reshape(d1velIndWW2B,(np.size(d1velIndWW2B),1))*d1nZ)
            st1bound[t]['Gamma'] = np.linalg.solve(d2aIJ, d1RHS)
            
            d2Gamma = st1bound[t]['Gamma'].reshape([iElemJ,iElemI],order='F')
            d2GammaCor = np.hstack([d2Gamma[:,0].reshape([iElemJ,1]), d2Gamma[:,1:] - d2Gamma[:,0:-1]])
            st1bound[t]['GammaCor'] = d2GammaCor.reshape([1,np.size(d2Gamma)],order='F')
            
            #lift
            d2dX = d2VRgridX[1:,0:-1]-d2VRgridX[0:-1,0:-1]
            d2dY = d2VRgridY[1:,0:-1]-d2VRgridY[0:-1,0:-1]
            d2dZ = d2VRgridZ[1:,0:-1]-d2VRgridZ[0:-1,0:-1]
            d2eNorm = np.sqrt(np.square(d2dX)+np.square(d2dY)+np.square(d2dZ))
            d2e1 = d2dX/d2eNorm
            d2e2 = d2dY/d2eNorm
            d2e3 = d2dZ/d2eNorm
            
            d2VtotU = d1Q[0] + np.sum(d2VelIndUB2B,axis=1).reshape([iElemJ,iElemI],order='F') + d1velIndUW2B.reshape([iElemJ,iElemI],order='F')
            d2VtotV = d1Q[1] + np.sum(d2VelIndVB2B,axis=1).reshape([iElemJ,iElemI],order='F') + d1velIndVW2B.reshape([iElemJ,iElemI],order='F')
            d2VtotW = d1Q[2] + np.sum(d2VelIndWB2B,axis=1).reshape([iElemJ,iElemI],order='F') + d1velIndWW2B.reshape([iElemJ,iElemI],order='F')
                  
#            d2Fx = drho*(d2VtotV*d2e3 - d2VtotW*d2e2)*d2GammaCor
#            d2Fy = drho*(d2VtotW*d2e1 - d2VtotU*d2e3)*d2GammaCor
#            d2Fz = drho*(d2VtotU*d2e2 - d2VtotV*d2e1)*d2GammaCor
            
            d2Fx = drho*(d1Q[1]*d2e3 - d1Q[2]*d2e2)*d2GammaCor
            d2Fy = drho*(d1Q[2]*d2e1 - d1Q[0]*d2e3)*d2GammaCor
            d2Fz = drho*(d1Q[0]*d2e2 - d1Q[1]*d2e1)*d2GammaCor
           
#            st1bound[t]['L'] = d2Fy/np.cos(np.arctan(d2VtotV/d2VtotU))
#            st1bound[t]['L'] = d2Fy/np.cos(np.arctan(d1Q[1]/d1Q[0]))
            st1bound[t]['L'] = d2Fx*(-np.sin(np.deg2rad(d1alpha[i])))+d2Fy*np.cos(np.deg2rad(d1alpha[i]))
            st1bound[t]['Ltot'] = np.sum(np.sum(st1bound[t]['L'],axis=1)*np.diff(d1spanPos))
            st1bound[t]['Fax'] = np.sum(np.sum(d2Fx,axis=1)*np.diff(d1spanPos))
            
            #lift coefficient
            d1dc = np.interp(d1spanPos[0:-1]+0.5*np.diff(d1spanPos),d1zLE,d1chord)
            
#            st1bound[t]['Cl'] = np.sum(st1bound[t]['L']/(0.5*drho*(np.square(d2VtotU)+np.square(d2VtotV))),axis=1)/d1dc
            st1bound[t]['Cl'] = np.sum(st1bound[t]['L']/(0.5*drho*np.square(dQ)),axis=1)/d1dc            
                            
            #Drag
            d2DWU = np.reshape(np.sum(d2velDWU,axis=1),(iElemJ,iElemI),order='F')
            d2DWV = np.reshape(np.sum(d2velDWV,axis=1),(iElemJ,iElemI),order='F')
            d2DWW = np.reshape(np.sum(d2velDWW,axis=1),(iElemJ,iElemI),order='F')
            
            d2Dx = -drho*(d2DWV*d2e3 - d2DWW*d2e2)*d2GammaCor
            d2Dy = -drho*(d2DWW*d2e1 - d2DWU*d2e3)*d2GammaCor
            d2Dz = -drho*(d2DWU*d2e2 - d2DWV*d2e1)*d2GammaCor
            
            st1bound[t]['D'] = d2Dx/np.cos(np.arctan(d2VtotV/d2VtotU))
            st1bound[t]['Cd'] = np.sum(st1bound[t]['D']/(0.5*drho*(np.square(d2VtotU)+np.square(d2VtotV))),axis=1)/d1dc
        
#            print(str("{:.2f}".format(round(d1alpha[i],2)))+'\t'+
#                  str("{:.2f}".format(iiter))+'\t'+
#                  str("{:.8f}".format(round(st1bound[t]['Cl'][int(np.ceil(iElemJ/2))],8)))+'\t'+
#                  str("{:.8f}".format(round(st1bound[t]['Cd'][int(np.ceil(iElemJ/2))],8)))+'\t'+
#                  str("{:.8f}".format(np.mean([dres1,dres2,dres3]))))
            print(str("{:.2f}".format(round(d1alpha[i],2)))+'\t'+
                  str("{:.2f}".format(iiter))+'\t'+
                  str("{:.8f}".format(round(st1bound[t]['Cl'][int(np.floor(iElemJ/2))],8)))+'\t'+
                  str("{:.8f}".format(round(st1bound[t]['Cd'][int(np.floor(iElemJ/2))],8)))+'\t'+
                  str("{:.8f}".format(np.mean([dres1,dres2,dres3]))))


        if iiter > 3:
#            dres1 = (st1bound[t]['Cl'][int(np.ceil(iElemJ/2))]-
#                     st1bound[t-1]['Cl'][int(np.ceil(iElemJ/2))])/st1bound[t-1]['Cl'][int(np.ceil(iElemJ/2))]
#            dres2 = (st1bound[t]['Cl'][int(np.ceil(iElemJ/2))]-
#                     st1bound[t-2]['Cl'][int(np.ceil(iElemJ/2))])/st1bound[t-2]['Cl'][int(np.ceil(iElemJ/2))]
#            dres3 = (st1bound[t]['Cl'][int(np.ceil(iElemJ/2))]-
#                     st1bound[t-3]['Cl'][int(np.ceil(iElemJ/2))])/st1bound[t-3]['Cl'][int(np.ceil(iElemJ/2))]
            dres1 = (st1bound[t]['Cl'][int(np.floor(iElemJ/2))]-
                     st1bound[t-1]['Cl'][int(np.floor(iElemJ/2))])/st1bound[t-1]['Cl'][int(np.floor(iElemJ/2))]
            dres2 = (st1bound[t]['Cl'][int(np.floor(iElemJ/2))]-
                     st1bound[t-2]['Cl'][int(np.floor(iElemJ/2))])/st1bound[t-2]['Cl'][int(np.floor(iElemJ/2))]
            dres3 = (st1bound[t]['Cl'][int(np.floor(iElemJ/2))]-
                     st1bound[t-3]['Cl'][int(np.floor(iElemJ/2))])/st1bound[t-3]['Cl'][int(np.floor(iElemJ/2))]
        
        iiter = iiter + 1
        t = t + 1
    
    if iElemJ == 1:
        d1cl[i] = st1bound[t-1]['Cl'][int(np.floor(iElemJ/2))]
        d1cd[i] = st1bound[t-1]['Cd'][int(np.floor(iElemJ/2))]
        d1Gamma[i] = np.sum(np.reshape(st1bound[t-1]['GammaCor'],(iElemJ,iElemI),order='F'),axis=1)[int(np.floor(iElemJ/2))]
    else:  
        d1cl[i] = st1bound[t-1]['Cl'][int(np.ceil(iElemJ/2))]
        d1cd[i] = st1bound[t-1]['Cd'][int(np.ceil(iElemJ/2))]
        d1Gamma[i] = np.sum(np.reshape(st1bound[t-1]['GammaCor'],(iElemJ,iElemI),order='F'),axis=1)[int(np.ceil(iElemJ/2))]
    
#%% plotting
d2XFOIL = np.fromfile('NACA0012polarsInviscidXFOIL.txt',dtype=float,sep=' ')

#circulation
fig, axG = plt.subplots(1,1)
#axL.clear()
axG.plot(d1alpha, np.pi*np.sin(np.deg2rad(d1alpha))*dchord*dQ,label='Thin airfoil theory')
axG.plot(d1alpha, -d1Gamma,label='Lifting surface code')
axG.set_title('NACA 0012')
axG.set_xlabel('Angle of attack [deg]')
axG.set_ylabel('$\Gamma$ [$m^2/s$]')
axG.grid('both')
axG.legend()

fig = plt.figure()
ax = fig.add_subplot()
ax.plot(d1alpha,-d1Gamma/(np.pi*np.sin(np.deg2rad(d1alpha))*dchord*dQ))
#ax.plot(d1alpha,np.cos(np.deg2rad(d1alpha)))
ax.set_xlabel('Angle of attack [deg]')
ax.set_ylabel('$\Gamma_{LSC}/\Gamma_{analytical}$ [-]')
ax.set_title('NACA 0012')
ax.grid('both')


#lift coefficient
fig, axL = plt.subplots(1,1)
#axL.clear()
axL.plot(d1alpha, 2*np.pi*np.sin(np.deg2rad(d1alpha))*np.cos(np.deg2rad(dsweep)),label='Thin airfoil theory')
axL.plot(d1alpha, d1cl,label='Lifting surface code')
axL.plot(d2XFOIL[::7],d2XFOIL[1::7],label='XFOIL')
axL.set_title('NACA 0012')
axL.set_xlabel('Angle of attack [deg]')
axL.set_ylabel('$c_l$ [-]')
axL.grid('both')
axL.legend()

#fig = plt.figure()
#ax = fig.add_subplot()
#ax.plot(d1alpha,d2XFOIL[1::7]/d1cl)
#ax.set_xlabel('Angle of attack [deg]')
#ax.set_ylabel('$c_{l,XFOIL}/c_{l,LSC}$ [-]')
#ax.set_title('NACA 4412')
#ax.grid('both')

#axD.clear()
#axD.plot(d2XFOIL[::7],d2XFOIL[3::7])
#axD.plot(d1alpha, d1cd)
#axD.set_xlabel('Angle of attack [deg]')
#axD.set_ylabel('$c_d$ [-]')
#axD.grid('both')


#%%writing output
d2out = np.vstack([d1alpha,d1cl,d1cd]).T
#np.savetxt('NACA4412polarsRe1500000.txt',d2out)
