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

import modelbase


class PETC2014(modelbase.Model):
    def __init__(self, p, r):
        super().__init__(p)
        
        compounds = [
                "PQ",  # oxidised plastoquinone
                "PC",  # oxidised plastocyan
                "Fd",  # oxidised ferrodoxin
                "ATP",  # stromal concentration of ATP
                "NADPH",  # stromal concentration of NADPH
                "H",  # lumenal protons
                "LHC",  # non-phosphorylated antenna
                "Fluo",
                "Light"]
        
        self.set_cpds(compounds)

        # Algebraic modules
        self.add_algebraicModule(r.pqmoiety, "pq_alm", ["PQ"],  ["PQred"])
        self.add_algebraicModule(r.pcmoiety, "pc_alm", ["PC"],  ["PCred"])
        self.add_algebraicModule(r.fdmoiety, "fd_alm", ["Fd"],  ["Fdred"])
        self.add_algebraicModule(r.adpmoiety, "adp_alm", ["ATP"],  ["ADP"])
        self.add_algebraicModule(r.nadpmoiety, "nadp_alm", ["NADPH"], ["NADP"])

        # add light driven electron transport chain reaction rates        

        self.set_ratev('vPS2', r.vPS2, "PQ", "PQred", "LHC")
        self.set_ratev("vPS1", r.vPS1, "PC", "PCred", "Fd", "Fdred", "LHC")
        self.set_ratev("vPTOX", r.vPTOX, "PQ", "PQred") #oxygen as a function of time
        self.set_rate("vB6f", r.vB6f, "PQ", "PQred", "PC", "PCred", "H")
        self.set_ratev("vNDH", r.vNDH, "PQ")
        self.set_rate("vCyc", r.vCyc, "PQ", "Fdred")
        self.set_rate("vFNR", r.vFNR, "Fd", "Fdred", "NADPH", "NADP")
        self.set_rate("vLeak", r.vLeak, "H")
        self.set_rate("vSt12", r.vSt12, "PQ", "LHC")
        self.set_rate("vSt21", r.vSt21, "LHC")
        self.set_rate("vATPsynthase", r.vATPsynthase, "ATP", "ADP", "H")
        self.set_rate("vATPconsumption", r.vATPconsumption, "ATP")
        self.set_rate("vNADPHconsumption", r.vNADPHconsumption, "NADPH")
        self.set_ratev("fluorescence", r.fluorescence, "PQ", "PQred", "LHC")
        self.set_ratev("light", r.light)

        self.set_stoichiometry_byCpd("PQ", {"vPS2": -1, "vB6f":1, "vCyc": -1, "vPTOX": 1, "vNDH": -1})
        self.set_stoichiometry_byCpd("PC", {"vB6f": -2, "vPS1": 1})
        self.set_stoichiometry_byCpd("Fd", {"vPS1": -1, "vFNR": 2, "vCyc": 2})
        self.set_stoichiometry_byCpd("ATP", {"vATPsynthase": 1, "vATPconsumption": -1})
        self.set_stoichiometry_byCpd("NADPH", {"vFNR": 1, "vNADPHconsumption": -1})
        self.set_stoichiometry_byCpd("H", {"vPS2": 2/p.bH, "vB6f": 4/p.bH, "vATPsynthase": -p.HPR/p.bH, "vLeak": -1/p.bH})
        self.set_stoichiometry_byCpd("LHC", {"vSt12": -1, "vSt21": 1})


        # finally add fluorescence as a dynamic variable to monitor teh dynamics of photosynthesis
        self.set_stoichiometry_byCpd("Fluo", {"fluorescence": 1})
        self.set_stoichiometry_byCpd("Light", {"light":1})
