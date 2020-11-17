# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 16:43:20 2020

@author: fritzek
"""

import numpy as np

def vrGrid(d2gridX,d2gridY,d2gridZ):
    
    d2gridXRoll = np.roll(d2gridX,-1,axis=1)
    d2gridYRoll = np.roll(d2gridY,-1,axis=1)
    d2gridZRoll = np.roll(d2gridZ,-1,axis=1)
    
    #vortex rings
    d2VRgridX       = d2gridX       + 0.25*(d2gridXRoll - d2gridX)
    d2VRgridX[:,-1] = d2gridX[:,-1] + 0.25*(d2gridX[:,-1] - d2gridX[:,-2])
    d2VRgridY       = d2gridY       + 0.25*(d2gridYRoll - d2gridY)
    d2VRgridY[:,-1] = d2gridY[:,-1] + 0.25*(d2gridY[:,-1] - d2gridY[:,-2])
    d2VRgridZ       = d2gridZ       + 0.25*(d2gridZRoll - d2gridZ)
    d2VRgridZ[:,-1] = d2gridZ[:,-1] + 0.25*(d2gridZ[:,-1] - d2gridZ[:,-2])
        
    return d2VRgridX, d2VRgridY, d2VRgridZ

def nVec(d2gridX,d2gridY,d2gridZ):
    
    d2gridAX = d2gridX[0:-1,1:]
    d2gridAY = d2gridY[0:-1,1:]
    d2gridAZ = d2gridZ[0:-1,1:]
    
    d2gridBX = d2gridX[0:-1,0:-1]
    d2gridBY = d2gridY[0:-1,0:-1]
    d2gridBZ = d2gridZ[0:-1,0:-1]
    
    d2gridCX = d2gridX[1:,0:-1]
    d2gridCY = d2gridY[1:,0:-1]
    d2gridCZ = d2gridZ[1:,0:-1]
    
    d1gridDX = d2gridX[1:,1:]
    d1gridDY = d2gridY[1:,1:]
    d1gridDZ = d2gridZ[1:,1:]
    
    d2nX = (d2gridCY-d2gridAY)*(d1gridDZ-d2gridBZ)-(d2gridCZ-d2gridAZ)*(d1gridDY-d2gridBY)
    d2nY = (d2gridCZ-d2gridAZ)*(d1gridDX-d2gridBX)-(d2gridCX-d2gridAX)*(d1gridDZ-d2gridBZ)
    d2nZ = (d2gridCX-d2gridAX)*(d1gridDY-d2gridBY)-(d2gridCY-d2gridAY)*(d1gridDX-d2gridBX)
    
    return d2nX,d2nY,d2nZ

def cpGrid(d2gridX,d2gridY,d2gridZ):
    
    d2gridXRollX = np.roll(d2gridX,-1,axis=1)
    d2gridYRollX = np.roll(d2gridY,-1,axis=1)
    d2gridZRollX = np.roll(d2gridZ,-1,axis=1)
    
#    d1gridXRollX = np.reshape(d2gridXRollX,(1,np.size(d2gridXRollX)),order='F')
#    d1gridYRollX = np.reshape(d2gridYRollX,(1,np.size(d2gridYRollX)),order='F')
#    d1gridZRollX = np.reshape(d2gridZRollX,(1,np.size(d2gridZRollX)),order='F')
    
    d2gridXRollZ = np.roll(d2gridX,-1,axis=0)
    d2gridYRollZ = np.roll(d2gridY,-1,axis=0)
    d2gridZRollZ = np.roll(d2gridZ,-1,axis=0)
    
#    d1gridXRollZ = np.reshape(d2gridXRollZ,(1,np.size(d2gridXRollZ)),order='F')
#    d1gridYRollZ = np.reshape(d2gridYRollZ,(1,np.size(d2gridYRollZ)),order='F')
#    d1gridZRollZ = np.reshape(d2gridZRollZ,(1,np.size(d2gridZRollZ)),order='F')
    
    d2gridXRollXZ = np.roll(d2gridXRollX,-1,axis=0)
    d2gridYRollXZ = np.roll(d2gridYRollX,-1,axis=0)
    d2gridZRollXZ = np.roll(d2gridZRollX,-1,axis=0)
    
#    d1gridXRollXZ = np.reshape(d2gridXRollXZ,(1,np.size(d2gridXRollXZ)),order='F')
#    d1gridYRollXZ = np.reshape(d2gridYRollXZ,(1,np.size(d2gridYRollXZ)),order='F')
#    d1gridZRollXZ = np.reshape(d2gridZRollXZ,(1,np.size(d2gridZRollXZ)),order='F')
    
    d2CPgridX = 0.5*(d2gridX + 0.75*(d2gridXRollX - d2gridX) +
                     d2gridXRollZ + 0.75*(d2gridXRollXZ-d2gridXRollZ))
    d2CPgridX = d2CPgridX[0:-1,0:-1]
    
    d2CPgridY = 0.5*(d2gridY + 0.75*(d2gridYRollX - d2gridY) +
                     d2gridYRollZ + 0.75*(d2gridYRollXZ-d2gridYRollZ))
    d2CPgridY = d2CPgridY[0:-1,0:-1]
    
    d2CPgridZ = 0.5*(d2gridZ + 0.75*(d2gridZRollX - d2gridZ) +
                     d2gridZRollZ + 0.75*(d2gridZRollXZ-d2gridZRollZ))
    d2CPgridZ = d2CPgridZ[0:-1,0:-1]
    
#    d1CPgridX = np.reshape(d2CPgridX,(1,np.size(d2CPgridX)),order='F')
#    d1CPgridY = np.reshape(d2CPgridY,(1,np.size(d2CPgridY)),order='F')
#    d1CPgridZ = np.reshape(d2CPgridZ,(1,np.size(d2CPgridZ)),order='F')
    
    return d2CPgridX,d2CPgridY,d2CPgridZ

