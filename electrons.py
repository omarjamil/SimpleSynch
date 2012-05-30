# This file is part of SimpleSynch.
#
# SimpleSynch is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SimpleSynch is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
#
# Author: Omar Jamil

import math

class Electrons:
    def __init__(self, g_min, g_max, p_law, ene_dens):
        #energy density in GeV/m^3
        self.value = self.kappa(g_min, g_max, p_law, ene_dens)
        
    def kappa(self, g_min, g_max, p_law, ene_dens):
        mc2 = 8.187111168006824e-14
        e_density = (ene_dens)/mc2

        if e_density == 0.0:
            kappa = 0.0
        elif g_max == 1.0:
            kappa = 0.0
        else:
            if p_law == 2.0:
                kappa = e_density/((math.log(g_max) - math.log(g_min)) + ((1./(g_max))-(1./(g_min))))
            else:
                temp1 = 1./(((1./(2.-p_law))*((g_max**(2.-p_law)) -
                                       (g_min**(2.-p_law)))) -
                   ((1./(1.-p_law))*((g_max**(1.-p_law)) - 
                                       (g_min**(1.-p_law)))))
                kappa = (mc2)**(p_law-1.) * e_density * temp1
        return kappa
