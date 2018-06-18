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

from dotmap import DotMap

class Parameters(DotMap):
    def __init__(self):
        super().__init__()  
        # pool sizes
        self.PSIItot = 2.5 # [mmol/molChl] total concentration of PSII
        self.PSItot = 2.5
        self.PQtot = 17.5 # [mmol/molChl]
        self.PCtot = 4. # Bohme1987 but other sources give different values - seems to depend greatly on organism and conditions
        self.Fdtot = 5. # Bohme1987
        self.Ctot = 2.5 #source unclear (Schoettler says 0.4...?, but plausible to assume that complexes (PSII,PSI,b6f) have approx. same abundance)
        self.NADPtot = 25. # estimate from ~ 0.8 mM, Heineke1991
        self.APtot = 60 # [mmol/molChl] Bionumbers ~2.55mM (=81mmol/molChl) (FIXME: Soma had 50)

        # parameters associated with photosystem II
        self.kH = 0.
        self.kH0 = 5e8 # base quenchingself. after calculation with Giovanni
        self.kF = 6.25e7 # 6.25e7 fluorescence 16ns
        self.k1 = 5e9 # excitation of Pheo / charge separation 200ps
        self.k1rev = 1.e10
        self.k2 = 5.e9 # original 5e9 (charge separation limiting step ~ 200ps) - made this faster for higher Fs fluorescence

        # parameters associated with photosystem I
        self.kStt7 = 0.0035 # [s-1] fitted to the FM dynamics
        self.kPph1 = 0.0013 # [s-1] fitted to the FM dynamics
        self.KM_ST = 0.2 # Switch point (half-activity of Stt7) for 20% PQ oxidised (80% reduced)
        self.n_ST = 2. # Hill coefficient of 4 -> 1/(2.5^4)~1/40 activity at PQox=PQred
        self.staticAntI = 0.2     # corresponds to PSI - LHCI supercomplex, when chlorophyll decreases more relative fixed antennae
        self.staticAntII = 0     # corresponds to PSII core


         # ATP and NADPH parameters
        self.kATPsynth = 20.    # taken from MATLAB
        self.kATPcons = 10.     # taken from MATLAB
        self.ATPcyt = 0.5       # only relative levels are relevant (normalised to 1) to set equilibrium
        self.Pi_mol = 0.01
        self.DeltaG0_ATP = 30.6 # 30.6kJ/mol / RT
        self.HPR = 14./3.  #Vollmar et al. 2009 (after Zhu et al. 2013)
        self.kNADPHcons = 15. # taken from MATLAB
        self.NADPHcyt = 0.5 # only relatice levels

        # pH and protons
        self.pHstroma = 7.8
        self.kLeak = 0.010 # [1/s] leakage rate -- inconsistency with Kathrine
        self.bH = 100. # proton buffer: ratio total / free protons

        # rate constants
        self.kPQred = 250. # [1/(s*(mmol/molChl))]
        self.kCytb6f = 2.5 # a rough estimate: transfer PQ->cytf should be ~10ms
        self.kPTOX = .01 # ~ 5 electrons / seconds. This gives a bit more (~20)
        self.kPCox = 2500. # a rough estimate: half life of PC->P700 should be ~0.2ms
        self.kFdred = 2.5e5 # a rough estimate: half life of PC->P700 should be ~2micro-s
        self.kcatFNR = 500. # Carrillo2003 (kcat~500 1/s)
        self.kcyc = 1.

        self.O2ext = 8. # corresponds to 250 microM corr. to 20%
        self.kNDH = .002 # re-introduce e- into PQ pool. Only positive for anaerobic (reducing) condition
        self.kNh = 0.05
        self.kNr = 0.004
        self.NPQsw = 5.8
        self.nH = 5.

        self.EFNR = 3. # Bohme1987
        self.KM_FNR_F = 1.56 # corresponds to 0.05 mM (Aliverti1990)
        self.KM_FNR_N = 0.22 # corresponds to 0.007 mM (Shin1971 Aliverti2004)

        # standard redox potentials (at pH=0) in V
        self.E0_QA = -0.140
        self.E0_PQ = 0.354
        self.E0_cytf = 0.350
        self.E0_PC = 0.380
        self.E0_P700 = 0.480
        self.E0_FA = -0.550
        self.E0_Fd = -0.430
        self.E0_NADP = -0.113

        # physical constants
        self.F = 96.485 # Faraday constant
        self.R = 8.3e-3 # universal gas constant
        self.T = 298. # Temperature in K - for now assumed to be constant at 25 C

        # light
        self.PFD = 0.
        self.Ton = 360
        self.Toff = 1800
        self.dT=120

        self.ox = True

