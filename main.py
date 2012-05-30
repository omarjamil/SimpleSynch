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

from photons import *
from electrons import *
from fluxDoppler import *
from param import *

if __name__  == "__main__":

#frequency values
    nu_val = [0] * grid_n
#iNu grid
    i_nu = [0] * grid_n

    electrons = Electrons(g_min, g_max, p, ene_dens)
    kappa = electrons.value
    photons = Photons(nu_min, nu_max, grid_n, nu_val, i_nu, kappa, B, p, radius)
#doppler and redshift the flux and frequencies
    fluxShift = FluxDoppler(BLF, theta, redshift, radius, nu_val, i_nu)
    
#output the results to a file
    f = open('output.dat', 'w')
    for i,j in map(None,nu_val,i_nu):
        s=str(i)+'\t'+str(j)+'\n'
        f.write(s)

    f.close()
    
    
