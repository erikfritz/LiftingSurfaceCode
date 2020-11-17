# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 16:00:56 2020

@author: fritzek
"""

import numpy as np

def vortexRing(d1P,d1A,d1B,d1C,d1D,dGamma):
        
    d1Vel1 = vortexLine(d1P,d1B,d1A,dGamma)
    d1Vel2 = vortexLine(d1P,d1C,d1B,dGamma)
    d1Vel3 = vortexLine(d1P,d1D,d1C,dGamma)
    d1Vel4 = vortexLine(d1P,d1A,d1D,dGamma)
    
    d1VelInd = d1Vel1+d1Vel2+d1Vel3+d1Vel4
    d1VelDW = d1Vel1+d1Vel3
    
    return d1VelInd,d1VelDW;
    
def vortexLine(d1P0,d1P1,d1P2,dGamma):
    
    d1r01 = d1P0-d1P1
    d1r02 = d1P0-d1P2
    
    d1n12 = np.cross(d1r01,d1r02)
    dnSquared = np.dot(d1n12,d1n12)
    
    dr01len = np.sqrt(np.sum(np.square(d1r01), axis=0))
    dr02len = np.sqrt(np.sum(np.square(d1r02), axis=0))
    
    if (dr01len<1e-10) or (dr02len<1e-10) or (dnSquared<1e-10):
        d1VelInd = np.array([0,0,0])
    else:
        d1r12 = d1P2-d1P1
        ddot1 = np.dot(d1r12,d1r01)
        ddot2 = np.dot(d1r12,d1r02)
        
        dcoeff = dGamma/(4*np.pi*dnSquared)*(ddot1/dr01len-ddot2/dr02len)
        
        d1VelInd = dcoeff*d1n12
        
    return d1VelInd

def vortexRingVec(d2evalPX,d2evalPY,d2evalPZ,d2inflPX,d2inflPY,d2inflPZ,d1Gamma):
        
    d1evalPX = np.reshape(d2evalPX,(1,np.size(d2evalPX)),order='F')
    d1evalPY = np.reshape(d2evalPY,(1,np.size(d2evalPY)),order='F')
    d1evalPZ = np.reshape(d2evalPZ,(1,np.size(d2evalPZ)),order='F')
        
    iElemI = np.size(d2inflPX,axis=1)-1
    iElemJ = np.size(d2inflPX,axis=0)-1
        
    d1gridAX = d2inflPX[0:-1,1:].reshape((1,iElemI*iElemJ),order='F')
    d1gridAY = d2inflPY[0:-1,1:].reshape((1,iElemI*iElemJ),order='F')
    d1gridAZ = d2inflPZ[0:-1,1:].reshape((1,iElemI*iElemJ),order='F')
    
    d1gridBX = d2inflPX[0:-1,0:-1].reshape((1,iElemI*iElemJ),order='F')
    d1gridBY = d2inflPY[0:-1,0:-1].reshape((1,iElemI*iElemJ),order='F')
    d1gridBZ = d2inflPZ[0:-1,0:-1].reshape((1,iElemI*iElemJ),order='F')
    
    d1gridCX = d2inflPX[1:,0:-1].reshape((1,iElemI*iElemJ),order='F')
    d1gridCY = d2inflPY[1:,0:-1].reshape((1,iElemI*iElemJ),order='F')
    d1gridCZ = d2inflPZ[1:,0:-1].reshape((1,iElemI*iElemJ),order='F')
    
    d1gridDX = d2inflPX[1:,1:].reshape((1,iElemI*iElemJ),order='F')
    d1gridDY = d2inflPY[1:,1:].reshape((1,iElemI*iElemJ),order='F')
    d1gridDZ = d2inflPZ[1:,1:].reshape((1,iElemI*iElemJ),order='F')
    
    d2Vel1U,d2Vel1V,d2Vel1W = vortexLineVec(d1evalPX,d1evalPY,d1evalPZ,
                                            d1gridBX,d1gridBY,d1gridBZ,
                                            d1gridAX,d1gridAY,d1gridAZ,
                                            d1Gamma)
    
    d2Vel2U,d2Vel2V,d2Vel2W = vortexLineVec(d1evalPX,d1evalPY,d1evalPZ,
                                            d1gridCX,d1gridCY,d1gridCZ,
                                            d1gridBX,d1gridBY,d1gridBZ,
                                            d1Gamma)
    
    d2Vel3U,d2Vel3V,d2Vel3W = vortexLineVec(d1evalPX,d1evalPY,d1evalPZ,
                                            d1gridDX,d1gridDY,d1gridDZ,
                                            d1gridCX,d1gridCY,d1gridCZ,
                                            d1Gamma)
    
    d2Vel4U,d2Vel4V,d2Vel4W = vortexLineVec(d1evalPX,d1evalPY,d1evalPZ,
                                            d1gridAX,d1gridAY,d1gridAZ,
                                            d1gridDX,d1gridDY,d1gridDZ,
                                            d1Gamma)
    
    d2VelIndU = d2Vel1U + d2Vel2U + d2Vel3U + d2Vel4U
    d2VelIndV = d2Vel1V + d2Vel2V + d2Vel3V + d2Vel4V
    d2VelIndW = d2Vel1W + d2Vel2W + d2Vel3W + d2Vel4W
    
    d2VelDWU = d2Vel1U + d2Vel3U
    d2VelDWV = d2Vel1V + d2Vel3V
    d2VelDWW = d2Vel1W + d2Vel3W

    return d2VelIndU,d2VelIndV,d2VelIndW,d2VelDWU,d2VelDWV,d2VelDWW

def vortexLineVec(d1P0X,d1P0Y,d1P0Z,d1P1X,d1P1Y,d1P1Z,d1P2X,d1P2Y,d1P2Z,d1Gamma):
    
    dcut = 1E-6
    
    d2r01X = d1P0X.T - d1P1X
    d2r01Y = d1P0Y.T - d1P1Y
    d2r01Z = d1P0Z.T - d1P1Z
    
    d2r02X = d1P0X.T - d1P2X
    d2r02Y = d1P0Y.T - d1P2Y
    d2r02Z = d1P0Z.T - d1P2Z
    
#    d1r12X = d1P2X - d1P1X
#    d1r12Y = d1P2Y - d1P1Y
#    d1r12Z = d1P2Z - d1P1Z
    
    d1r12X = d1P1X - d1P2X
    d1r12Y = d1P1Y - d1P2Y
    d1r12Z = d1P1Z - d1P2Z

    d2nX = d2r01Y*d2r02Z-d2r01Z*d2r02Y
    d2nY = d2r01Z*d2r02X-d2r01X*d2r02Z
    d2nZ = d2r01X*d2r02Y-d2r01Y*d2r02X
    
    d2n2 = np.square(d2nX)+np.square(d2nY)+np.square(d2nZ)
    d2r01norm = np.sqrt(np.square(d2r01X)+np.square(d2r01Y)+np.square(d2r01Z))
    d2r02norm = np.sqrt(np.square(d2r02X)+np.square(d2r02Y)+np.square(d2r02Z))
    
    b2n2 = d2n2 <= dcut
    b2r01norm = d2r01norm <= dcut
    b2r02norm = d2r02norm <= dcut
    
    b2tot = b2n2 | b2r01norm | b2r02norm
    
    d2dot1 = (np.tile(d1r12X,[np.size(d2r01X,axis=0),1])*d2r01X + 
              np.tile(d1r12Y,[np.size(d2r01Y,axis=0),1])*d2r01Y + 
              np.tile(d1r12Z,[np.size(d2r01Z,axis=0),1])*d2r01Z)
    d2dot2 = (np.tile(d1r12X,[np.size(d2r02X,axis=0),1])*d2r02X + 
              np.tile(d1r12Y,[np.size(d2r02Y,axis=0),1])*d2r02Y + 
              np.tile(d1r12Z,[np.size(d2r02Z,axis=0),1])*d2r02Z)
    
    d2Gamma = np.tile(d1Gamma.T,[np.size(d2r01X,axis=0),1])
    
    d2K = d2Gamma/(4*np.pi*d2n2)*(d2dot1/d2r01norm-d2dot2/d2r02norm)
    d2K = np.nan_to_num(d2K)
    d2K = 1*~b2tot*d2K
    
    d2Vel1U = d2nX*d2K
    d2Vel1V = d2nY*d2K
    d2Vel1W = d2nZ*d2K

    return d2Vel1U,d2Vel1V,d2Vel1W