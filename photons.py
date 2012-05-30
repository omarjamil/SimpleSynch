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
from synchrotron import *

class Photons:
    def __init__(self, nu_min, nu_max, grid_n, nu_val, i_nu, kappa, B, p, radius):
        length = radius
        self.freq_range(nu_min, nu_max, grid_n, nu_val)
        Synchrotron(kappa, p, B, nu_val, i_nu, length)

#logarithmic frequency range
    def freq_range(self, nu_min, nu_max, grid_n, nu_val):
        abscissa = nu_min
        r = (nu_max/nu_min)**(1./grid_n)
        r2 = math.sqrt(r)
        r -= 1.0
        for i in range(grid_n):
            nu_val[i] = r2*abscissa
            abscissa += r*abscissa
        

