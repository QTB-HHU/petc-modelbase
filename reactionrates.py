#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 08:20:58 2018

The mathematical model of the photosynthetic electron transport chain defines methods to calculate reaction rates
and set of ten differential equations based on the model published by Ebenhoeh et al. in 2014

Copyright (C) 2014-2018  Anna Matuszyńska, Oliver Ebenhöh

This program is free software: you can redistribute and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with this program (license.txt).  If not, see <http://www.gnu.org/licenses/>.
"""
__author__ = "Anna Matuszyńska"
__copyright__ = "Copyright 2018, Heinrich-Heine University Dusseldorf"
__credits__ = ["Anna Matuszynska", "Oliver Ebenhoeh"]
__maintainer__ = "Anna Matuszynska"
__email__ = "Anna.Matuszynska@uni-duesseldorf.de"
__status__ = "Production/Stable"

import numpy as np
import math
import warnings
# ================================================= #
# Photosynthetic Electron Transport Reaction Rates  #
# ================================================= #

class Reactions:

    def fluorescence(self, p, P, Pred, LHC, **kwargs):
        cs = self.crossSection(p, LHC)
        B = self.ps2states(p, P, Pred, LHC, **kwargs)
        fluo = cs * ((p.kF / (p.kF + p.kH0 + p.k2)) * B[0] + (p.kF / (p.kF + p.kH0)) * B[2])
        return fluo

    
    # ====================================================================== #
    # Composed parameters #
    # ====================================================================== #
    def RT(self,p):
        return p.R*p.T

    def dG_pH(self,p):
        return np.log(10)*self.RT(p)

    def Hstroma(self,p):
        return 3.2e4*10**(-p.pHstroma) 

    def kProtonation(self,p):
        return 4e-3 / self.Hstroma(p)

    def pH(self,hlf):
        if hlf<=0:
            warnings.warn("H became negative {} and was set to 1e-30".format(hlf))
            hlf=1e-30       
        return (-math.log(hlf*2.5e-4)/math.log(10))

    def pHinv(self,pH):
        return (4e3*10**(-pH))

    def Keq_PQred(self,p):
        DG1 = -p.E0_QA * p.F
        DG2 = -2 * p.E0_PQ * p.F
        DG = -2 * DG1 + DG2 + 2 * p.pHstroma * self.dG_pH(p)
        K = np.exp(-DG/self.RT(p))
        return K

    def Keq_cyc(self,p):
        DG1 = -p.E0_Fd * p.F
        DG2 = -2 * p.E0_PQ * p.F
        DG = -2 * DG1 + DG2 + 2 * self.dG_pH(p) * p.pHstroma
        K = np.exp(-DG/self.RT(p))
        return K

    def Keq_cytfPC(self,p):
        DG1 = -p.E0_cytf * p.F
        DG2 = -p.E0_PC * p.F
        DG = -DG1 + DG2
        K = np.exp(-DG/self.RT(p))
        return K

    def Keq_FAFd(self,p):
        DG1 = -p.E0_FA * p.F
        DG2 = -p.E0_Fd * p.F
        DG = -DG1 + DG2
        K = np.exp(-DG/self.RT(p))
        return K

    def Keq_PCP700(self,p):
        DG1 = -p.E0_PC * p.F
        DG2 = -p.E0_P700 * p.F
        DG = -DG1 + DG2
        K = np.exp(-DG/self.RT(p))
        return K

    def Keq_NDH(self,p):
        DG1 = -2 * p.E0_NADP * p.F
        DG2 = -2 * p.E0_PQ * p.F
        DG = -DG1 + DG2 + self.dG_pH(p) * p.pHstroma
        K = np.exp(-DG/self.RT(p))
        return K

    def Keq_FNR(self,p):
        DG1 = -p.E0_Fd * p.F
        DG2 = -2 * p.E0_NADP * p.F
        DG = -2 * DG1 + DG2 + self.dG_pH(p) * p.pHstroma
        K = np.exp(-DG/self.RT(p))
        return K

    def Keq_ATP(self,p, pH):
        DG = p.DeltaG0_ATP - self.dG_pH(p) * p.HPR * (p.pHstroma - pH)
        Keq = p.Pi_mol * np.exp(-DG/self.RT(p))
        return Keq

    def Keq_cytb6f(self,p, pH):
        DG1 = -2 * p.F * p.E0_PQ
        DG2 = -p.F * p.E0_PC
        DG = - (DG1 + 2*self.dG_pH(p) * pH) + 2 * DG2 + 2*self.dG_pH(p) * (p.pHstroma - pH)
        Keq = np.exp(-DG/self.RT(p))
        return Keq

    # ====================================================================== #
    # Conserved quantities -> for algebraic module #
    # ====================================================================== #
    def pqmoiety(self, p, PQ):
        return p.PQtot - PQ

    def pcmoiety(self, p, PC):
        return p.PCtot - PC

    def fdmoiety(self, p, Fd):
        return p.Fdtot - Fd

    def adpmoiety(self, p, ATP):
        return p.APtot - ATP

    def nadpmoiety(self, p, NADPH):
        return p.NADPtot - NADPH

    # ====================================================================== #
    # Dumb light function
    # ====================================================================== #

#    def light(self, p, Ton=270, Toff=900, Tfl=120, **kwargs):
#        '''
#        :return: light intensity at certain point of time. 
#        Typical PAM light function
#        '''
#        return 5000 * ((kwargs["t"] % Tfl) < 0.8) + \
#                p.PFD * (((kwargs["t"] > Ton) and (kwargs["t"] < Toff)) and \
#                         (kwargs["t"] % Tfl) >= 0.8)

    def light(self, p, **kwargs):
        '''
        :return: light intensity at certain point of time. 
        Typical PAM light function
        '''
        if kwargs["t"]%p.dT <=0.8:
            return 5000
        elif ((kwargs["t"] > p.Ton) and (kwargs["t"] < p.Toff) and \
              (kwargs["t"]%p.dT>0.8)):
            return p.PFD
        else:
            return 0.0000001


#    def light(self, p, **kwargs):
#        '''
#        :return: light intensity at certain point of time. 
#        Typical PAM light function
#        '''
#        return p.PFD
      

    # ====================================================================== #
    # Reaction rates
    # ====================================================================== #
    def ps2states(self, p, PQox, PQred, LHC, **kwargs):
        """ 
        QSSA, calculates the states of the photosystem II
        accepts values:
            Pox: oxidised fraction of the PQ pool (PQH2)
            Q: quencher
            L: light, int or array of the n x 1 dimension, that gives the light intensity

            returns:
            B: array of arrays with the states of PSII; rows: time, columns states: 1 and 3 excited
        """
        L = self.LII(p, LHC,**kwargs)

        Q = 0#self.vQuencher4states(p, Psbs, Viola)

        k2 = p.k2
        kF = p.kF
        kH = p.kH0 #+ p.kH * Q
        k3p = p.kPQred * PQox
        k3m = p.kPQred * PQred / self.Keq_PQred(p)
        M = np.array([[-L-k3m, kH+kF,       k3p, 0],
                      [L,      -(kH+kF+k2), 0,   0],
                      [0,      0,           L,   -(kH+kF)],
                      [1,      1,           1,   1]])

        A = np.array([0, 0, 0, p.PSIItot])
        B = np.linalg.solve(M, A)
        return B

    def ps1states(self, p, PC, PCred, Fd, Fdred, LHC,**kwargs):
        """ 
        QSSA calculates open state of PSI
        depends on reduction states of plastocyanin and ferredoxin
        C = [PC], F = [Fd] (ox. forms)
        accepts: light, y as an array of arrays
        returns: array of PSI open
        """
        L = self.LI(p, LHC,**kwargs)

        A1 = p.PSItot / (1 + L/(p.kFdred * Fd) + (1 + Fdred/(self.Keq_FAFd(p) * Fd))
                          * (PC/(self.Keq_PCP700(p) * PCred)
                             + L/(p.kPCox * PCred))
        )
        return A1

    ###############################################################################
    # method to calculate cross sections
    ###############################################################################

    def crossSection(self, p, LHC):
        """ calculates the cross section of PSII """
        cs = p.staticAntII + (1 - p.staticAntII - p.staticAntI) * LHC
        return cs

    def LII(self,p,LHC,**kwargs):
        return self.crossSection(p, LHC) * self.light(p, **kwargs)

    def LI(self,p,LHC,**kwargs):
        return (1-self.crossSection(p, LHC)) * self.light(p, **kwargs)

    ###############################################################################
    # Reaction rates
    ###############################################################################
     
    def vPS2(self, p, P, Pred, LHC, **kwargs):
        """ reaction rate constant for photochemistry """
        #Q = self.quencher(p,Q,H)
        B = self.ps2states(p, P, Pred, LHC, **kwargs)
        v = p.k2 * B[1] / 2
        return v

    def vPS1(self, p, PC, PCred, Fd, Fdred, LHC, **kwargs):
        """ reaction rate constant for open PSI """
        L = self.LI(p, LHC, **kwargs)
        A = self.ps1states(p, PC, PCred, Fd, Fdred, LHC, **kwargs)
        v = L * A
        return v
    
    def oxygen(self, p, **kwargs):
        """ return oxygen and NDH concentration as a function of time
        used to simulate anoxia conditions as in the paper"""
        if p.ox == True:
            ''' by default we assume constant oxygen supply'''
            return p.O2ext, p.kNDH
        else:
            if kwargs['t']<= p.Ton or kwargs['t']>=p.Toff:
                return p.O2ext, 0
            else:
                return 0, p.kNDH

    def vPTOX(self, p, Pox, Pred, **kwargs):
        """ calculates reaction rate of PTOX """
        v = Pred * p.kPTOX * self.oxygen(p, **kwargs)[0] 
        return v

    def vNDH(self, p, Pox, **kwargs):
        """ 
        calculates reaction rate of PQ reduction under absence of oxygen
        can be mediated by NADH reductase NDH
        """
        v = self.oxygen(p, **kwargs)[1] * Pox
        return v

    def vB6f(self, p, Pox, Pred, PC, PCred, H):
        """ calculates reaction rate of cytb6f """
        ph = self.pH(H)
        Keq = self.Keq_cytb6f(p,ph)
        v = max(p.kCytb6f * (Pred * PC**2 - (Pox * PCred**2)/Keq), -p.kCytb6f)
        return v

    def vCyc(self, p, Pox, Fdred):
        """
        calculates reaction rate of cyclic electron flow
        considered as practically irreversible
        """
        v = p.kcyc * ((Fdred**2) * Pox)
        return v

    def vFNR(self, p, Fd, Fdred, NADPH, NADP):
        """
        Reaction rate mediated by the Ferredoxin—NADP(+) reductase (FNR)
        Kinetic: convenience kinetics Liebermeister and Klipp, 2006
        Compartment: lumenal side of the thylakoid membrane
        Units:
        Reaction rate: mmol/mol Chl/s
        [F], [Fdred] in mmol/mol Chl/s
        [NADPH] in mM
        """
        fdred = Fdred/p.KM_FNR_F
        fdox = Fd/p.KM_FNR_F
        nadph = (NADPH)/p.KM_FNR_N  # NADPH requires conversion to mmol/mol of chlorophyll
        nadp = (NADP)/p.KM_FNR_N
        v = (p.EFNR * p.kcatFNR *
            ((fdred**2) * nadp - ((fdox**2) * nadph) / self.Keq_FNR(p)) /
            ((1+fdred+fdred**2) * (1+nadp) + (1+fdox+fdox**2) * (1+nadph) - 1))
        return v

    def vLeak(self, p, Hlf):
        """ 
        rate of leak of protons through the membrane
        """
        v = p.kLeak * (Hlf - self.pHinv(p.pHstroma))
        return v

    def vSt12(self,p, Pox,Ant):
        """ 
        reaction rate of state transitions from PSII to PSI
        Ant depending on module used corresponds to non-phosphorylated antennae
        or antennae associated with PSII
        """
        kKin = p.kStt7 * ( 1 / (1 + ((Pox /p.PQtot)/p.KM_ST)**p.n_ST))
        v = kKin * Ant
        return v

    def vSt21(self,p, Ant):
        """
        reaction rate of state transitions from PSI to PSII
        """
        v = p.kPph1 * (1 - Ant)
        return v

    def vATPsynthase(self, p, ATP, ADP, H):
        """
        Reaction rate of ATP production
        Kinetic: simple mass action with PH dependant equilibrium
        Compartment: lumenal side of the thylakoid membrane
        """
        ph = self.pH(H)
        v = p.kATPsynth * (ADP - ATP / self.Keq_ATP(p, ph))
        return v
        
        
    def vATPconsumption(self, p, A):
        v = p.kATPcons * A
        return v


    def vNADPHconsumption(self, p, N):
        v = p.kNADPHcons * N
        return v


    def vQuencherBase(self, p, Q, H):
        vNPQdyn = H**p.nH / (H ** p.nH + self.pHinv(p.NPQsw)**p.nH)
        v = (1-Q) * p.kNh * vNPQdyn - p.kNr * Q
        return v    