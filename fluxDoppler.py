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
# along with SimpleSynch.  If not, see <http://www.gnu.org/licenses/>.
#
# Author: Omar Jamil

import math

class FluxDoppler:
    def __init__(self, BLF, theta, redshift, radius, nu_val, i_nu):
        DL_Mpc = self.nedWright(redshift)
        doppVal = self.bulkDopplerFactor(BLF, theta)
        freqShift = doppVal/(1.+redshift)
        area = 4. * math.pi * radius**2
        dF = self.dopplerFlux(doppVal, DL_Mpc, area)
        
#Shift the frequency and flux according to the source ditance, speed
        for i in range(len(nu_val)):
            nu_val[i] *= freqShift
            i_nu[i] *= dF 

#source flux
    def dopplerFlux(self, doppVal, DL_Mpc, area):
        sourceDist = DL_Mpc * 3.08568025e22;
        denom = 4. * math.pi * (sourceDist**2);
        val = ((area)/denom)*(doppVal**3);
        return val;
      
#source doppler factor
    def bulkDopplerFactor(self, lorentzF, viewTheta):
        beta = math.sqrt(1. - (1./(lorentzF)**2));
        dopplerF = 1./(lorentzF * (1-beta*(math.cos(viewTheta))));
        return dopplerF;

#source luminosity distance
    def nedWright(self, redshift):
        # Reference: Wright, E. L., 2006, PASP, 118, 1711.
        z = redshift
        n=1000	# number of points in integrals
        c = 299792.458 # velocity of light in km/sec
        Tyr = 977.8 #coefficent for converting 1/H into Gyr
        h0 = 0.7
        H0 = h0*100.	# Hubble constant
        Omega_M = 0.270
        Omega_Lambda = 0.730
  
        # Densities
        WM = Omega_M		# Omega(matter)
        WV = Omega_Lambda	# Omega(vacuum) or lambda
        h = H0/100.	# H0/100
        WR = 4.165e-5/(h*h)	# Omega(radiation), includes 3 massless neutrino species, T0 = 2.72528
        WK = 1-WM-WR-WV	# Omega curvaturve = 1-Omega(total)
  
        a = 1.0	# the scale factor of the Universe
        az = 1.0/(1.+1.0*z)
  
        age = 0
        for i in range(n):
            a = az*(float(i)+0.5)/float(n)
            adot = math.sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
            age = age + 1/adot
        
        zage = az*age/n
        # correction for annihilations of particles not present now like e+/e-
        # added 13-Aug-03 based on T_vs_t.f
        lpz = math.log((1.+1.0*z))/math.log(10.0)
        dzage = 0.
        if (lpz >  7.500):
            dzage = 0.002 * (lpz -  7.500)
        if (lpz >  8.000):
            dzage = 0.014 * (lpz -  8.000) +  0.001
        if (lpz >  8.500):
            dzage = 0.040 * (lpz -  8.500) +  0.008
        if (lpz >  9.000):
            dzage = 0.020 * (lpz -  9.000) +  0.028
        if (lpz >  9.500):
            dzage = 0.019 * (lpz -  9.500) +  0.039
        if (lpz > 10.000):
            dzage = 0.048
        if (lpz > 10.775):
            dzage = 0.035 * (lpz - 10.775) +  0.048
        if (lpz > 11.851):
            dzage = 0.069 * (lpz - 11.851) +  0.086
        if (lpz > 12.258):
            dzage = 0.461 * (lpz - 12.258) +  0.114
        if (lpz > 12.382):
            dzage = 0.024 * (lpz - 12.382) +  0.171
        if (lpz > 13.055):
            dzage = 0.013 * (lpz - 13.055) +  0.188
        if (lpz > 14.081):
            dzage = 0.013 * (lpz - 14.081) +  0.201
        if (lpz > 15.107):
            dzage = 0.214
        zage = zage*(10.0**dzage)
        zage_Gyr = (Tyr/H0)*zage
  
        DTT = 0.0
        DCMR = 0.0
        # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
        for i in range(n):
            a = az+(1.-az)*(float(i)+0.5)/float(n)
            adot = math.sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
            DTT = DTT + 1./adot
            DCMR = DCMR + 1./(a*adot)
            
        DTT = (1-az)*DTT/float(n)
        DCMR = (1-az)*DCMR/float(n)
  
        age = DTT+zage
        age_Gyr = age*(Tyr/H0)
	
        DTT_Gyr = (Tyr/H0)*DTT
	
        DCMR_Gyr = (Tyr/H0)*DCMR
        DCMR_Mpc = (c/H0)*DCMR
  
        # tangential comoving distance
        ratio = 1.00
        x = math.sqrt(math.fabs(WK))*DCMR
        if x > 0.1: 
            if WK>0:
                ratio = 0.5*(math.exp(x)-math.exp(-x))/x
            else:
                ratio = math.sin(x)/x
                y = ratio*DCMR
               
        else:	
            y = x*x
        # statement below fixed 13-Aug-03 to correct sign error in expansion
            if WK < 0:
                y = -y
                ratio = 1. + y/6. + y*y/120.
                y= ratio*DCMR
                
        DCMT = y
        DA = az*DCMT
        DA_Mpc = (c/H0)*DA
        kpc_DA = DA_Mpc/206.264806
        DA_Gyr = (Tyr/H0)*DA
        DL = DA/(az*az)
        DL_Mpc = (c/H0)*DL
        DL_Gyr = (Tyr/H0)*DL
  
        # printf("Distance parameters: age_Gyr %e zage_Gyr %e DTT_Gyr %e DA_Mpc %e kpc_DA %e DL_Mpc %e \n", *age_Gyr,*zage_Gyr,*DTT_Gyr,*DA_Mpc,*kpc_DA,*DL_Mpc)
	
        return DL_Mpc



