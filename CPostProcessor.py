# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 08:12:43 2020

@author: fritzek
"""
#%%import packages
import numpy as np
import matplotlib.pyplot as plt 

#%%class definition
class PostProcessor:
    
    def plotWake(self,iElemI,iElemJ,turbine,tPlot):
        import matlab.engine
        eng = matlab.engine.start_matlab()
        
        d3bladeX = np.ndarray([iElemJ+1,iElemI+1,turbine.nBlades])
        d3bladeY = np.ndarray([iElemJ+1,iElemI+1,turbine.nBlades])
        d3bladeZ = np.ndarray([iElemJ+1,iElemI+1,turbine.nBlades])
        d3wakeX = np.ndarray([iElemJ+1,tPlot+1,turbine.nBlades])
        d3wakeY = np.ndarray([iElemJ+1,tPlot+1,turbine.nBlades])
        d3wakeZ = np.ndarray([iElemJ+1,tPlot+1,turbine.nBlades])
        d3wakeG = np.ndarray([iElemJ,tPlot,turbine.nBlades])
        
        for i in range(0,turbine.nBlades):
            d3bladeX[:,:,i] = turbine.blade[i].d2X[tPlot]
            d3bladeY[:,:,i] = turbine.blade[i].d2Y[tPlot]
            d3bladeZ[:,:,i] = turbine.blade[i].d2Z[tPlot]
            d3wakeX[:,:,i] = turbine.wake[i].d2X[tPlot]
            d3wakeY[:,:,i] = turbine.wake[i].d2Y[tPlot]
            d3wakeZ[:,:,i] = turbine.wake[i].d2Z[tPlot]
            
            d2Gamma = np.reshape(turbine.wake[i].Gamma[tPlot],[iElemJ,tPlot],order='F')
            d3wakeG[:,:,i] = d2Gamma
            
        R = eng.plotGeometry(matlab.double(d3bladeX.tolist()),
                             matlab.double(d3bladeY.tolist()),
                             matlab.double(d3bladeZ.tolist()),
                             matlab.double(d3wakeX.tolist()),
                             matlab.double(d3wakeY.tolist()),
                             matlab.double(d3wakeZ.tolist()),
                             matlab.double(d3wakeG.tolist()),
                             turbine.hubHeight)
        
        return eng

    def plotCirculation(self,iElemI,iElemJ,turbine,tPlot):
        fig = plt.figure('Circulation')
        ax = fig.add_subplot()
        ax.clear()
        for i in range(0,turbine.nBlades):
            d1span = turbine.blade[i].d1spanPos[0:-1] + 0.5*np.diff(turbine.blade[i].d1spanPos) - turbine.blade[i].bladeRoot
            d1Gamma = np.sum(np.reshape(turbine.blade[i].GammaCor[tPlot],(iElemJ,iElemI),order='F'),axis=1)
            lbl = 'Blade %i' % i
            ax.plot(d1span, d1Gamma, label= lbl)
        ax.set_xlabel('Span [m]')
        ax.set_ylabel('Circulation $\Gamma$ [$\mathrm{m}^2/\mathrm{s}$]')
        ax.grid(b=True,which='major')
        ax.grid(b=True,which='minor')
        ax.legend()
        
    def plotLift(self,turbine,tPlot):
        fig = plt.figure('Lift')
        ax = fig.add_subplot()
        ax.clear()
        for i in range(0,turbine.nBlades):
            d1span = turbine.blade[i].d1spanPos[0:-1] + 0.5*np.diff(turbine.blade[i].d1spanPos) - turbine.blade[i].bladeRoot
            lbl = 'Blade %i' % i
            ax.plot(d1span, turbine.blade[i].d1L[tPlot], label= lbl)
        ax.set_xlabel('Span [m]')
        ax.set_ylabel('Lift $L$ [$\mathrm{N}/\mathrm{m}$]')
        ax.grid(which='both')
        ax.legend()
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
        
    def plotCl(self,turbine,tPlot):
        fig = plt.figure('Lift coefficient')
        ax = fig.add_subplot()
        ax.clear()
        for i in range(0,turbine.nBlades):
            d1span = turbine.blade[i].d1spanPos[0:-1] + 0.5*np.diff(turbine.blade[i].d1spanPos) - turbine.blade[i].bladeRoot
            lbl = 'Blade %i' % i
            ax.plot(d1span, turbine.blade[i].d1cl[tPlot], label= lbl)
        ax.set_xlabel('Span [m]')
        ax.set_ylabel('Lift coefficient $c_l$ [-]')
        ax.grid(which='both')
        ax.legend()
        
    def plotPitchMoment(self,turbine,tPlot):
        fig = plt.figure('Pitching Moment')
        ax = fig.add_subplot()
        ax.clear()
        for i in range(0,turbine.nBlades):
            dazPitchAx  = np.deg2rad(turbine.bladePos[tPlot][i])
#            d2azForce   = np.pi - np.arctan(-turbine.blade[i].d2Ycp[tPlot]/(turbine.blade[i].d2Zcp[tPlot]-turbine.hubHeight-turbine.zNac2hub))
            d2azForce   = 2*np.arctan(-turbine.blade[i].d2Ycp[tPlot]/(np.sqrt(np.square(turbine.blade[i].d2Ycp[tPlot])+np.square(turbine.blade[i].d2Zcp[tPlot]-turbine.hubHeight-turbine.zNac2hub))+turbine.blade[i].d2Zcp[tPlot]-turbine.hubHeight-turbine.zNac2hub))
            d2rForce    = np.sqrt(np.square(turbine.blade[i].d2Ycp[tPlot])+np.square(turbine.blade[i].d2Zcp[tPlot]-turbine.hubHeight-turbine.zNac2hub))
            d2rPitchAx  = d2rForce*np.cos(d2azForce-dazPitchAx)
            
            d2rPitchX = d2rPitchAx*0
            d2rPitchY = -d2rPitchAx*np.sin(dazPitchAx)
            d2rPitchZ = d2rPitchAx*np.cos(dazPitchAx)
            
            d2lX = turbine.blade[i].d2Xcp[tPlot] - turbine.xNac2hub - d2rPitchX
            d2lY = turbine.blade[i].d2Ycp[tPlot] - d2rPitchY
            d2lZ = turbine.blade[i].d2Zcp[tPlot] - turbine.hubHeight - turbine.zNac2hub - d2rPitchZ
            
            d2MX = turbine.blade[i].d2Fy[tPlot]*d2lZ - turbine.blade[i].d2Fz[tPlot]*d2lY
            d2MY = turbine.blade[i].d2Fz[tPlot]*d2lX - turbine.blade[i].d2Fx[tPlot]*d2lZ
            d2MZ = turbine.blade[i].d2Fx[tPlot]*d2lY - turbine.blade[i].d2Fy[tPlot]*d2lX
            
            d2Mpitch = - d2MY*np.sin(dazPitchAx) + d2MZ*np.cos(dazPitchAx)
            d1Mpitch = np.sum(d2Mpitch,axis=1)
            dMpitch = np.sum(d1Mpitch)
            
            d1span = turbine.blade[i].d1spanPos[0:-1] + 0.5*np.diff(turbine.blade[i].d1spanPos) - turbine.blade[i].bladeRoot
            lbl = 'Blade %i' % i
            ax.plot(d1span,d1Mpitch,marker='.',label=lbl)

        ax.set_title('Total pitching moment $M_{p,tot}=%f$ Nm' % dMpitch)
        ax.set_xlabel('Span [m]')
        ax.set_ylabel('Pitching moment $M_{pitch}$ [Nm]')
        ax.grid(which='both')
        ax.legend()
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
        
        return d1Mpitch
        
    def plotGlobalLoads(self,sim,turbine,tPlot):
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.clear()
        d1torque = np.zeros_like(sim.t)
        d1P = np.zeros_like(sim.t)
        for i in range(1,sim.tsteps):
            d1torque[i] = turbine.dtorque[i]
            d1P[i] = turbine.dP[i]
        ax.plot(sim.t, d1torque,'k')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Torque $T$ [Nm]')
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
        ax.grid(which='both')
        ax2 = ax.twinx()
        ax2.plot(sim.t, d1P,'--k')
        ax2.set_ylabel('Power $P$ [W]')
        ax2.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
        
    def plotVelPlane(self,sim,turbine,tPlot,dX,dY,dZ,sDir,dRes):
        from vortexRing import vortexRingVec
        
        minX = minY = minZ =  1000.0
        maxX = maxY = maxZ = -1000.0
        for i in range(0,turbine.nBlades):
            minX = np.min([minX,np.min(turbine.blade[i].d2X[tPlot]),np.min(turbine.wake[i].d2X[tPlot])])
            minY = np.min([minY,np.min(turbine.blade[i].d2Y[tPlot]),np.min(turbine.wake[i].d2Y[tPlot])])
            minZ = np.min([minZ,np.min(turbine.blade[i].d2Z[tPlot]),np.min(turbine.wake[i].d2Z[tPlot])])
            maxX = np.max([maxX,np.max(turbine.blade[i].d2X[tPlot]),np.max(turbine.wake[i].d2X[tPlot])])
            maxY = np.max([maxY,np.max(turbine.blade[i].d2Y[tPlot]),np.max(turbine.wake[i].d2Y[tPlot])])
            maxZ = np.max([maxZ,np.max(turbine.blade[i].d2Z[tPlot]),np.max(turbine.wake[i].d2Z[tPlot])])
            
        if sDir == 'X':
            iDir = 0
            iN1 = int(np.round((maxY-minY)/dRes,0))
            dresY = (maxY-minY)/iN1
            iN2 = int(np.round((maxZ-minZ)/dRes,0))
            dresZ = (maxZ-minZ)/iN2
            
            d2gridX = np.tile(dX,(iN2+1,iN1+1))
            d2gridY = np.tile(np.arange(minY,maxY+0.01,dresY),(iN2+1,1))
            d2gridZ = np.tile(np.arange(minZ,maxZ+0.01,dresZ).reshape(iN2+1,1),(1,iN1+1))
            
        elif sDir == 'Y':
            iDir = 1
            iN1 = int(np.round((maxX-minX)/dRes,0))
            dresX = (maxX-minX)/iN1
            iN2 = int(np.round((maxZ-minZ)/dRes,0))
            dresZ = (maxZ-minZ)/iN2
            
            d2gridX = np.tile(np.arange(minX,maxX+0.01,dresX),(iN2+1,1))
            d2gridY = np.tile(dY,(iN2+1,iN1+1))
            d2gridZ = np.tile(np.arange(minZ,maxZ+0.01,dresZ).reshape(iN2+1,1),(1,iN1+1))
            
        elif sDir == 'Z':
            iDir = 2
            iN1 = int(np.round((maxX-minX)/dRes,0))
            dresX = (maxX-minX)/iN1
            iN2 = int(np.round((maxY-minY)/dRes,0))
            dresY = (maxY-minY)/iN2
            
            d2gridX = np.tile(np.arange(minX,maxX+0.01,dresX),(iN2+1,1))
            d2gridY = np.tile(np.arange(minY,maxY+0.01,dresY).reshape(iN2+1,1),(1,iN1+1))
            d2gridZ = np.tile(dZ,(iN2+1,iN1+1))
            
        else:
            print('Wrong directional input.')
        
        #B2P
        for i in range(0,turbine.nBlades):
            tmpU,tmpV,tmpW,_,_,_ = vortexRingVec(d2gridX,d2gridY,d2gridZ,
                                                 turbine.blade[i].d2X[tPlot],turbine.blade[i].d2Y[tPlot],turbine.blade[i].d2Z[tPlot],
                                                 turbine.blade[i].Gamma[tPlot])
            if i == 0:
                d2velIndB2Pu = tmpU
                d2velIndB2Pv = tmpV
                d2velIndB2Pw = tmpW
            else:
                d2velIndB2Pu = d2velIndB2Pu + tmpU
                d2velIndB2Pv = d2velIndB2Pv + tmpV
                d2velIndB2Pw = d2velIndB2Pw + tmpW

            tmpU,tmpV,tmpW,_,_,_ = vortexRingVec(d2gridX,d2gridY,d2gridZ,
                                                 turbine.wake[i].d2X[tPlot],turbine.wake[i].d2Y[tPlot],turbine.wake[i].d2Z[tPlot],
                                                 turbine.wake[i].Gamma[tPlot])
            if i == 0:
                d2velIndW2Pu = tmpU
                d2velIndW2Pv = tmpV
                d2velIndW2Pw = tmpW
            else:
                d2velIndW2Pu = d2velIndW2Pu + tmpU
                d2velIndW2Pv = d2velIndW2Pv + tmpV
                d2velIndW2Pw = d2velIndW2Pw + tmpW

        d2U = (sim.d2wind[0,1] + np.sum(d2velIndB2Pu,axis=1) + np.sum(d2velIndW2Pu,axis=1)).reshape((iN1+1,iN2+1)).T
        d2V = (sim.d2wind[0,2] + np.sum(d2velIndB2Pv,axis=1) + np.sum(d2velIndW2Pv,axis=1)).reshape((iN1+1,iN2+1)).T
        d2W = (sim.d2wind[0,3] + np.sum(d2velIndB2Pw,axis=1) + np.sum(d2velIndW2Pw,axis=1)).reshape((iN1+1,iN2+1)).T
        
        sXlabel = ['y [m]','x [m]','x [m]']
        sYlabel = ['z [m]','z [m]','y [m]']
        d1loc = [dX,dY,dZ]
        l1grid = [d2gridX,d2gridY,d2gridZ]
        l1plot = [[1,2],[0,2],[0,1]]
        
        fig, axs = plt.subplots(2, 2)
        
        contU = axs[0,0].contourf(l1grid[l1plot[iDir][0]],l1grid[l1plot[iDir][1]],np.clip(d2U,np.mean(d2U)-4*np.std(d2U),np.mean(d2U)+4*np.std(d2U)),20)
        axs[0,0].set_xlabel(sXlabel[iDir])
        axs[0,0].set_ylabel(sYlabel[iDir])
        axs[0,0].set_title('U at '+sDir+' = '+str(d1loc[iDir])+'m')
        axs[0,0].set_aspect('equal', 'box')
        cbarU = fig.colorbar(contU,ax=axs[0,0])
        contU.set_clim(np.mean(d2U)-4*np.std(d2U),np.mean(d2U)+4*np.std(d2U))
        cbarU.set_label('U [m/s]')

        contV = axs[0,1].contourf(l1grid[l1plot[iDir][0]],l1grid[l1plot[iDir][1]],np.clip(d2V,np.mean(d2V)-4*np.std(d2V),np.mean(d2V)+4*np.std(d2V)),20)
        axs[0,1].set_xlabel(sXlabel[iDir])
        axs[0,1].set_ylabel(sYlabel[iDir])
        axs[0,1].set_title('V at '+sDir+' = '+str(d1loc[iDir])+'m')
        axs[0,1].set_aspect('equal', 'box')
        cbarV = fig.colorbar(contV,ax=axs[0,1])
        cbarV.set_label('V [m/s]')

        contW = axs[1,0].contourf(l1grid[l1plot[iDir][0]],l1grid[l1plot[iDir][1]],np.clip(d2W,np.mean(d2W)-4*np.std(d2W),np.mean(d2W)+4*np.std(d2W)),20)
        axs[1,0].set_xlabel(sXlabel[iDir])
        axs[1,0].set_ylabel(sYlabel[iDir])
        axs[1,0].set_title('W at '+sDir+' = '+str(d1loc[iDir])+'m')
        axs[1,0].set_aspect('equal', 'box')
        cbarW = fig.colorbar(contW,ax=axs[1,0])
        cbarW.set_label('W [m/s]')
        
        l1vel = [d2U,d2V,d2W]
        l1circ = [[2,1],[0,2],[1,0]]
        l1diff = [[1,0],[0,1],[1,0]]
        
        d2dvel1 = np.diff(l1vel[l1circ[iDir][0]],axis=l1diff[iDir][0])
        d2dvel2 = np.diff(l1vel[l1circ[iDir][1]],axis=l1diff[iDir][1])
        
        d2dx1 = np.diff(l1grid[l1circ[iDir][1]],axis=l1diff[iDir][0])
        d2dx2 = np.diff(l1grid[l1circ[iDir][0]],axis=l1diff[iDir][1])
        
        if l1diff[iDir][0]:
            d2dvel1 = np.vstack([d2dvel1[:,0],d2dvel1.T]).T
            d2dvel2 = np.vstack([d2dvel2[0,:],d2dvel2])
            
            d2dx1 = np.vstack([d2dx1[:,0],d2dx1.T]).T
            d2dx2 = np.vstack([d2dx2[0,:],d2dx2])
        else:
            d2dvel1 = np.vstack([d2dvel1[0,:],d2dvel1])
            d2dvel2 = np.vstack([d2dvel2[:,0],d2dvel2.T]).T
            
            d2dx1 = np.vstack([d2dx1[0,:],d2dx1])
            d2dx2 = np.vstack([d2dx2[:,0],d2dx2.T]).T
    
        d2circ = d2dvel1/d2dx1 - d2dvel2/d2dx2
        
        contC = axs[1,1].contourf(l1grid[l1plot[iDir][0]],l1grid[l1plot[iDir][1]],np.clip(d2circ,np.mean(d2circ)-4*np.std(d2circ),np.mean(d2circ)+4*np.std(d2circ)),20)
        axs[1,1].set_xlabel(sXlabel[iDir])
        axs[1,1].set_ylabel(sYlabel[iDir])
        axs[1,1].set_title('$\omega_'+sDir+'$ at '+sDir+' = '+str(d1loc[iDir])+'m')
        axs[1,1].set_aspect('equal', 'box')
        cbarC = fig.colorbar(contC,ax=axs[1,1])
        cbarC.set_label('$\omega$ [$s^-1$]')
        
        fig.tight_layout()

