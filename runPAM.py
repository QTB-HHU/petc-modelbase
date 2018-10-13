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
# Import built-in libraries
import instantiate
import numpy as np

#For graphical display
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()


def pamsimulation(s, y0, Tmax, Ton, Toff, photon=100, dT=120, Tflash=0.8):  
    '''actual PAM simulation with dark adaptation 
    by default: light-dark-light transition'''
    
    s.set_initial_value(y0, t0=0)
    t = 0
    y = y0
    
    s.model.par.update({'continuous':False, 'Ton':Ton, 'Toff': Toff, 'dT':dT, 'PFD':photon})
    
    s.set_initial_value(y0, t0=0)
    try:
        '''if Assimulo present'''
        s.integrator.atol = 1.e-8
    except:
        pass
    t = 0
    y = y0

    while s.successful() and t < Tmax:             
        #turn on the saturating pulse of light of Tflash length
        if t%s.model.par.dT == 0:
            tfl = np.linspace(t, t+Tflash, 100)
            s.timeCourse(tfl, y)
        else:
            #switch on the light except for the dark period
            #t+dT-Tflash is the time to the next flash
            if t<= s.model.par.Ton or t>=s.model.par.Toff:
                tIfl = np.linspace(t, t+s.model.par.dT-Tflash, 10000)
                s.timeCourse(tIfl, y)
            else:
                #put the actinic light
                tIfl = np.linspace(t, t+s.model.par.dT-Tflash, 1000)
                s.timeCourse(tIfl, y)
        t = s.getT()[-1]
        y = s.getVarsByName(s.model.cpdNames)[-1]
        
    return s, s.model

if __name__=='__main__':
    p, rat, m, s = instantiate.instantiate()
    
    # Pass initial values
    init = {"PQ":m.par.PQtot/2,
            "PC":m.par.PCtot/2,
            "Fd":m.par.Fdtot/2,
            "ATP":m.par.APtot/3,
            "NADPH":m.par.NADPtot/2,
            "H (lumen)":rat.pHinv(7.2),
            "LHC":0.9,
            'Fluo':0.2,
            "Light":m.par.PFD}
    
    y0 = np.array(list(init.values()))   
    
    s, m = pamsimulation(s, y0, 1800, 270, 900, 100, 90, 0.8)


    plt.figure()
    fl = s.getRate('fluorescence')
    plt.plot(s.getT(), fl/max(fl))
    plt.xlabel('Time [s]')
    plt.ylabel('Fluorescence ~PSII [normalized]')

    plt.show()
