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

def anoxia(s, y0, Tmax=3600, Tflash=0.8):
    """Returns the solution of time integration with anoxia induced state transitions
    Anoxia is simulated by shutting down import reaction ssuch as PTOX and 
    setting the atmospheric oxygen to 0.
    
    Keyword arguments:
        Tmax -- time of experiment [in seconds]
        Ton -- time when the oxygen is taken [s]
        Toff -- time when the oxygen is put back into the system [s]
    """
    t = 0
    y = y0

    s.model.par.update({'PFD': 0, 'ox':False})
    try:
        '''only executed if assimulo is used'''
        s.integrator.atol = 1.e-8
    except:
        pass

    # =========================================== #
    # Anoxia experiment simulation in three parts #
    while s.successful() and t < Tmax:    
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
    init = {"PQ":m.par.PQtot,
            "PC":0.0202,
            "Fd":5.000,
            "ATP":0.0,
            "NADPH":0.0,
            "H (lumen)":rat.pHinv(7.2),
            "LHC":0.9,
            'Fluo':0,
            'Light':p.PFD}
        
    y0 = np.array(list(init.values()))   
       
    s, m = anoxia(s, y0)

    plt.figure()
    fl = s.getRate('fluorescence')
    plt.plot(s.getT(), fl/max(fl))
    plt.xlabel('Time [s]')
    plt.ylabel('Fluorescence ~PSII [normalized]')

    plt.show()