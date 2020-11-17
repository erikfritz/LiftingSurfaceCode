# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 14:50:38 2020

@author: fritzek
"""
#%% TO DO
"""
include:
    yaw angle
    tilt angle
    cone angle
"""

#%% import packages
import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#from drawnow import drawnow, figure
#from tabulate import tabulate
#from scipy import interpolate
from scipy.interpolate import griddata
import geomParser
from CBlade import Blade
from CWake import Wake
#import time
#%% class definition
class Turbine:
    def __init__(self, dcinput):
        self.cone           = float(dcinput['CONEANGLE'])
        self.groundLvl      = float(dcinput['GROUNDLEVEL'])
        self.hubHeight      = float(dcinput['HUBHEIGHT'])
        self.nBlades        = int(dcinput['NROFBLADES'])
        self.rpm            = float(dcinput['RPM'])
        self.startAz        = float(dcinput['STARTAZIMUTH'])
        self.tilt           = float(dcinput['TILTANGLE'])
        self.towerBaseR     = float(dcinput['TOWERBASERADIUS'])
        self.towerTopR      = float(dcinput['TOWERTOPRADIUS'])
        self.xNac2hub       = float(dcinput['XNAC2HUB']) 
        self.yaw            = float(dcinput['YAWANGLE'])
        self.zNac2hub       = float(dcinput['ZNAC2HUB']) 
        self.dFax           = {}
        self.dtorque        = {}
        self.dP             = {}
        self.blade          = [None]*self.nBlades
        self.wake           = [None]*self.nBlades
        self.bladePos       = {}
        self.bladePos[0]    = (np.arange(0.,360.,360/self.nBlades)+self.startAz)%360
        
    def initBlades(self,dcinput,iElemI,iElemJ):
        for i in range(0,self.nBlades):
            self.blade[i] = Blade(self,i,dcinput)
            self.blade[i].iElemI = iElemI
            self.blade[i].iElemJ = iElemJ
            
    def discrBlade(self,sfolder,sairfoilFile,d2aeroprops,iElemI,iElemJ):
        #load airfoils and calculate camber surface
        fairfoils = open(sfolder+'/'+sairfoilFile,'r')
        l1airfoils = fairfoils.readlines()
                
        d1z = np.array(l1airfoils[0].replace('63.7hr_r','').split(),dtype='float') #TODO: soft code
        
        d2airfoils = np.empty([0,len(d1z)+1])
        for i in range(1,len(l1airfoils)):
            d1tmp = np.array(l1airfoils[i].split(),dtype='float') 
            d2airfoils = np.vstack([d2airfoils,d1tmp])
            
        iLE = np.argmin(d2airfoils[:,0])
        d2SS = d2airfoils[0:iLE+1,:]
        d2SS = d2SS[d2SS[:,0].argsort()]
        d2PS = d2airfoils[iLE:,:]
        
        d2SSinterp = np.empty([100,np.size(d2airfoils,axis=1)])
        d2PSinterp = np.empty([100,np.size(d2airfoils,axis=1)])
        
        d2SSinterp[:,0] = np.linspace(0,1,100)
        d2PSinterp[:,0] = np.linspace(0,1,100)
        
        for i in range(1,np.size(d2airfoils,axis=1)):
            d2SSinterp[:,i] = np.interp(np.linspace(0,1,100),d2SS[:,0],d2SS[:,i])
            d2PSinterp[:,i] = np.interp(np.linspace(0,1,100),d2PS[:,0],d2PS[:,i])
            
        d2camb = 0.5*(d2SSinterp[:,1:]+d2PSinterp[:,1:])
        
        #calculate panels
        d1x     = d2aeroprops[:,5]
        d1y     = d2aeroprops[:,6]
        d1c     = d2aeroprops[:,1]
        d1t     = np.deg2rad(-d2aeroprops[:,3])
        d1c14   = d2aeroprops[:,4]
        
        d1LEx = d1x + d1c*(0.5 - (25 + d1c14)/100)*np.sin(d1t)
        d1LEy = d1y - d1c*(0.5 - (25 + d1c14)/100)*np.cos(d1t)
        
        d2z = np.tile(d1z,(100,1))
        d2y = np.tile(np.linspace(0,1,100),(len(d1z),1)).T
        d2gridP = np.vstack([d2z.reshape(1,-1),d2y.reshape(1,-1)]).T
        d2gridV = d2camb.reshape(1,-1).T
        
        d1zPan = min(d1z)+(max(d1z)-min(d1z))*(np.cos(np.linspace(-np.pi,0,iElemJ+1))+1)/2
        d2zPan = np.tile(d1zPan,(iElemI+1,1))
        d2yInterp = np.tile(np.linspace(0,1,iElemI+1),(iElemJ+1,1)).T
        
        d2gridNew = np.vstack([d2zPan.reshape(1,-1),d2yInterp.reshape(1,-1)]).T
        d1cambInterp = griddata(d2gridP,d2gridV,d2gridNew,method='linear')
        d2cambInterp = d1cambInterp.reshape(iElemI+1,iElemJ+1)
        
        d1LExInterp = np.interp(d1zPan,d1z,d1LEx)
        d1LEyInterp = np.interp(d1zPan,d1z,d1LEy)
        d1cInterp = np.interp(d1zPan,d1z,d1c)
        d1tInterp = np.interp(d1zPan,d1z,d1t)
        
        d2xPan = np.empty([iElemI+1,iElemJ+1])
        d2yPan = np.empty([iElemI+1,iElemJ+1])
        for i in range(0,np.size(d2cambInterp,axis=1)):
            d2cambTmp = np.vstack([np.linspace(0,1,iElemI+1),d2cambInterp[:,i]]).T
            d2rot = [[np.cos(-d1tInterp[i]),np.sin(-d1tInterp[i])],[-np.sin(-d1tInterp[i]),np.cos(-d1tInterp[i])]]
            d2cambRot = np.matmul(d2cambTmp,d2rot)
            d2xPan[:,i] = d1LExInterp[i]+d2cambRot[:,1]*d1cInterp[i]
            d2yPan[:,i] = d1LEyInterp[i]+d2cambRot[:,0]*d1cInterp[i]
            
        d2xVR,d2yVR,d2zVR = geomParser.vrGrid(d2xPan.T,d2yPan.T,d2zPan.T)
        d2xCP,d2yCP,d2zCP = geomParser.cpGrid(d2xPan.T,d2yPan.T,d2zPan.T)
        
        for i in range(0,self.nBlades):
            dpitch = np.deg2rad(self.blade[i].pitch)
            d2xPitch = np.cos(dpitch)*d2xVR-np.sin(dpitch)*d2yVR
            d2yPitch = np.sin(dpitch)*d2xVR+np.cos(dpitch)*d2yVR
            d2zPitch = d2zVR
            d2xCPpitch = np.cos(dpitch)*d2xCP-np.sin(dpitch)*d2yCP
            d2yCPpitch = np.sin(dpitch)*d2xCP+np.cos(dpitch)*d2yCP
            d2zCPpitch = d2zCP
            
            daz = np.deg2rad(self.bladePos[0][i])
            d2xPitchAz = d2xPitch
            d2yPitchAz = np.cos(daz)*d2yPitch-np.sin(daz)*d2zPitch
            d2zPitchAz = np.sin(daz)*d2yPitch+np.cos(daz)*d2zPitch
            d2xCPpitchAz = d2xCPpitch
            d2yCPpitchAz = np.cos(daz)*d2yCPpitch-np.sin(daz)*d2zCPpitch
            d2zCPpitchAz = np.sin(daz)*d2yCPpitch+np.cos(daz)*d2zCPpitch
            
            d1translate = np.array([self.xNac2hub,0,self.hubHeight+self.zNac2hub])
            d2xPitchAzTransl = d2xPitchAz+d1translate[0]
            d2yPitchAzTransl = d2yPitchAz
            d2zPitchAzTransl = d2zPitchAz+d1translate[2]
            d2xCPpitchAzTransl = d2xCPpitchAz+d1translate[0]
            d2yCPpitchAzTransl = d2yCPpitchAz
            d2zCPpitchAzTransl = d2zCPpitchAz+d1translate[2]

            self.blade[i].d2X[0] = d2xPitchAzTransl
            self.blade[i].d2Y[0] = d2yPitchAzTransl
            self.blade[i].d2Z[0] = d2zPitchAzTransl
            self.blade[i].d2Xcp[0] = d2xCPpitchAzTransl
            self.blade[i].d2Ycp[0] = d2yCPpitchAzTransl
            self.blade[i].d2Zcp[0] = d2zCPpitchAzTransl
            self.blade[i].d1spanPos = d1zPan

    def initWake(self):
        for i in range(0,self.nBlades):
            self.wake[i] = Wake(self,i)
            
    def rotate(self,t,dtStep):
        for i in range(0,self.nBlades):
            d1translate = np.array([self.xNac2hub,0,self.hubHeight+self.zNac2hub])
            d2tmpX = self.blade[i].d2X[t-1] - d1translate[0]
            d2tmpY = self.blade[i].d2Y[t-1]
            d2tmpZ = self.blade[i].d2Z[t-1] - d1translate[2]
            d2tmpXcp = self.blade[i].d2Xcp[t-1] - d1translate[0]
            d2tmpYcp = self.blade[i].d2Ycp[t-1]
            d2tmpZcp = self.blade[i].d2Zcp[t-1] - d1translate[2]
            
            daz = np.deg2rad(self.rpm*6*dtStep)
            self.blade[i].d2X[t] = d2tmpX + d1translate[0]
            self.blade[i].d2Y[t] = np.cos(daz)*d2tmpY - np.sin(daz)*d2tmpZ
            self.blade[i].d2Z[t] = np.sin(daz)*d2tmpY + np.cos(daz)*d2tmpZ + d1translate[2]
            self.blade[i].d2Xcp[t] = d2tmpXcp + d1translate[0]
            self.blade[i].d2Ycp[t] = np.cos(daz)*d2tmpYcp - np.sin(daz)*d2tmpZcp
            self.blade[i].d2Zcp[t] = np.sin(daz)*d2tmpYcp + np.cos(daz)*d2tmpZcp + d1translate[2]
            
        self.bladePos[t] = (self.bladePos[t-1] + self.rpm*6*dtStep)%360
        
    def calcTurbineLoads(self,t):
        dAz = np.deg2rad(self.bladePos[t])
        for i in range(0,self.nBlades):
            d2Faz = - np.sin(dAz[i])*self.blade[i].d2Fz[t] - np.cos(dAz[i])*self.blade[i].d2Fy[t]
            d2Torque = d2Faz*np.sqrt(np.square(self.blade[i].d2Ycp[t])+np.square(self.blade[i].d2Zcp[t]-self.hubHeight-self.zNac2hub))
            dFax = self.blade[i].dFaxTot[t]
            
            if i == 0:
#                self.d1Faz[t] = np.sum(d2Faz,axis=1)
                self.dtorque[t] = np.sum(np.sum(d2Torque,axis=1))
                self.dFax[t] = dFax
            else:
#                self.d1Faz[t] = self.d1Faz[t] + np.sum(d2Faz,axis=1)
                self.dtorque[t] = self.dtorque[t] + np.sum(np.sum(d2Torque,axis=1))
                self.dFax[t] = self.dFax[t] + dFax
                
            self.dP[t] = self.dtorque[t]*2*np.pi/60*self.rpm