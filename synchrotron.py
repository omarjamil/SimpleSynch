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

#Longair (High energy astrophysics, vol 2)
#synchrotron treatment for power-law electrons followed

class Synchrotron:
    def __init__(self, kappa, p_index, B, nu_val, i_nu, length):
        a_p = self.ap(p_index)
        b_p = self.bp(p_index)
        for i in range(len(nu_val)):
            nu = nu_val[i]
            em = self.emiss(a_p, B, p_index, kappa, nu)
            ab = self.absor(kappa, B, p_index, b_p, nu)
            i_nu[i] += (em/(4.*math.pi*ab))*(1.- (math.exp(-1.*ab*length)))
          
#Synchrotron emissivity
    def emiss(self, a_p, B, p_index, kappa, nu):
        emiss_val = 2.344e-25 * a_p * B**((p_index+1.)/2.)*kappa*((1.253e37)/(nu))**((p_index-1.)/2.)
        return emiss_val

#Synchrotron absorption
    def absor(self, kappa, B, p_index, b_p, nu):
        absor_val = 3.354e-9*kappa*(B**((p_index+2.)*0.5))*((3.54e18)**(p_index))*b_p*nu**((p_index+4.)*(-0.5))
        return absor_val

    def ap(self, p_index):
        gam_num = math.gamma((p_index/4.) + (19./12.))*math.gamma(p_index/4. - (1./12.))*math.gamma(p_index/4.+(5./4.))
        gam_den = (p_index+1.)*math.gamma(p_index/4.+(7./4.))
        gam = gam_num/gam_den
        a_p_val = 0.5 * math.sqrt(math.pi) * gam
        return a_p_val

    def bp(self, p_index):
        gam_num = math.gamma((3.*p_index)*(1./12.))*math.gamma((3.*p_index)*(1./12.))*math.gamma((p_index*6.)*(0.25))
        gam_den = math.gamma((p_index+8.)*(0.25))
        gam = gam_num/gam_den
        b_p_val = 1./8.*(math.sqrt(math.pi)) * gam
        return b_p_val
