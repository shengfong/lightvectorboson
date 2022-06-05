# math
from math import pi, exp, log, log10, sqrt
# numpy
import numpy as np
# scipy
from scipy.linalg import expm
from scipy.interpolate import interp1d
from scipy.special import zeta
# abc
from abc import ABC, abstractmethod

# input
from acropolis.input import InputInterface, locate_sm_file
# params
from acropolis.params import zeta3
from acropolis.params import hbar, c_si, me, me2, mmu, mmu2, mpic, mpic2, mpi, mpi2, alpha, tau_t

# model
from acropolis.models import AbstractModel

import sys

class DecayVectorModel(AbstractModel):

    def __init__(self, mphi, tau, temp0, n0a, bree, brmumu, brpipi, brpia, br3pi):
        # Initialize the Input_Interface
        self._sII   = InputInterface( locate_sm_file() )

        # The mass of the decaying particle
        self._sMphi = mphi            # in MeV
        # The lifetime of the decaying particle
        self._sTau  = tau             # in s
        # The injection energy
        self._sE0   = self._sMphi/2.  # in MeV #Energy of electron/positron from V -> e+e-
        ###### second "delta" source with energy E1 < E0
        if mphi > mpi:
            self._sE1   = self._sMphi/2.*(1.-mpi2/self._sMphi**2.) #Energy of photon from V -> pi a
        else:
            self._sE1   = self._sE0

        # The number density of the mediator
        # (relative to photons) ...
        self._sN0a  = n0a
        # ... at T = temp0 ...
        self._sT0   = temp0           # in MeV
        # ... corresponding to t = t(temp0)
        self._st0   = self._sII.time(self._sT0)

        # The branching ratio into electron-positron pairs
        self._sBRee = bree
        # The branching ratio into muon-antimuon pairs
        self._sBRmumu = brmumu
        # The branching ratio into pi+ pi- pairs
        self._sBRpipi = brpipi
        # The branching ratio into pi0 photon pairs
        self._sBRpia = brpia
        # The branching ratio into pi0 pi+ p-
        self._sBR3pi = br3pi

        # INTERPOLATION OF PHOTON AND ELECTRON/POSITION SPECTRA ################
        self._filename_gamma = []
        self._filename_e = []
        self._num = 100 #Number of points in log10 scale of energy of gamma/electron/positron for each MV
        self._filename_gamma.append("spec_data/spec_mumu_g.dat") # The gamma spectrum (excluding FSR) for V -> mu+ mu-: (MV,E,dN_gamma/dE)
        self._filename_gamma.append("spec_data/spec_pipi_g.dat") # The gamma spectrum (excluding FSR) for V -> pi+ pi-: (MV,E,dN_gamma/dE)
        self._filename_gamma.append("spec_data/spec_3pi_g.dat") # The gamma spectrum for V -> pi0 pi+ pi-: (MV,E,dN_gamma/dE)
        self._filename_e.append("spec_data/spec_mumu_e.dat") # The electron/positon spectrum for V -> mu+ mu-: (MV,E,dN_e/dE)
        self._filename_e.append("spec_data/spec_pipi_e.dat") # The electron/positon spectrum for V -> pi+ pi-: (MV,E,dN_e/dE)
        self._filename_e.append("spec_data/spec_3pi_e.dat") # The electron/positon spectrum for V -> pi0 pi+ pi-: (MV,E,dN_e/dE)
        
        #Read the corresponding index of MV 
        #Important: Use the same grid of MV for both gamma and electron/positron spectra
        MV = self._sMphi
        iV = []

        #print(mphi,bree,brmumu,brpipi,brpia,br3pi)
        self._MVspec = []
        for i in range(3):
            #skip the first row for labels
            data = np.loadtxt(self._filename_gamma[i], skiprows=1, usecols = (0)) 
            MV_i = []
            for j in range(int(data.shape[0]/self._num)):
                k = self._num*j
                MV_i.append(data[k])
            iV.append(np.argmin( np.abs( MV - np.array(MV_i) ))*self._num) #Identity the index of closest mV
            #It is possible that the MV in the spectrum is slightly larger or smaller than the input MV
            #If larger, no problem
            #If smaller, it might try to read the spectrum out of range of energy
            #So we set MV always equal to that of the spectrum 
            self._MVspec.append(MV_i[int(iV[i]/self._num)])
            #print(MV,MV_i[int(iV[i]/self._num)])
           
        #print(self._MVspec)
            
        
        self._spec_gamma = []
        self._spec_e = []
        for i in range(3):
            #skip the first row for labels
            data_tmp = np.loadtxt(self._filename_gamma[i], skiprows=1 + iV[i], usecols = (1,2), max_rows = self._num)
            data  = np.where(data_tmp > 1.0e-20, data_tmp, 1.0e-20) #If smaller then 1.e-20, set to 1.e-20
            EE_read = np.log10(data[:,0]) # [MeV]
            dNdE_read = np.log10(data[:,1]) # [-]
            self._spec_gamma.append(interp1d(EE_read, dNdE_read, kind='linear'))
            
            #skip the first row for labels
            data_tmp = np.loadtxt(self._filename_e[i], skiprows=1 + iV[i], usecols = (1,2), max_rows = self._num)
            data  = np.where(data_tmp > 1.0e-20, data_tmp, 1.0e-20) #If smaller then 1.e-20, set to 1.e-20
            EE_read = np.log10(data[:,0]) # [MeV]
            dNdE_read = np.log10(data[:,1]) # [-]
            self._spec_e.append(interp1d(EE_read, dNdE_read, kind='linear'))
            
            

        # Call the super constructor
        super(DecayVectorModel, self).__init__(self._sE0, self._sE1, self._sII)


    # DEPENDENT QUANTITIES ##############################################################

    def _number_density(self, T):
        sf_ratio = self._sII.scale_factor(self._sT0)/self._sII.scale_factor(T)

        delta_t = self._sII.time(T) - self._st0
        n_gamma = (2.*zeta3)*(self._sT0**3.)/(pi**2.)
    
        return self._sN0a * n_gamma * sf_ratio**3. * exp( -delta_t/self._sTau )


    # ABSTRACT METHODS ##################################################################

    def _temperature_range(self):
        # The number of degrees-of-freedom to span
        mag = 2.
        # Calculate the approximate decay temperature
        Td = self._sII.temperature( self._sTau )
        # Calculate Tmin and Tmax from Td
        Td_ofm = log10(Td)
        # Here we choose -1.5 (+0.5) orders of magnitude
        # below (above) the approx. decay temperature,
        # since the main part happens after t = \tau
        gstar = 10.633 #relativistics degree of freedom at T=1 MeV
        conv = np.pi**4/(45*zeta(3))*gstar #conversion factor s(T)/ngamma(T)
        if ( self._sTau < 1.e9 and self._sN0a*self._sMphi/conv < 1.e-9 ): 
            Tmin = 10.**(Td_ofm - 3.*mag/4.)
            Tmax = 10.**(Td_ofm + 1.*mag/4.)
        else:    
            Tmin = max(10.**(Td_ofm - 1.*mag),7.e-7) #Let us stop at the matter-radiation equality 7*10^{-7} MeV
            Tmax = 10.**(Td_ofm + 4.*mag/2.)
            if Tmax <= Tmin:
                Tmax = 10.**2.*Tmin
            
        return (Tmin, Tmax)


    def _source_photon_0(self, T):
        return 0.
        
    def _source_photon_1(self, T):
        return self._sBRpia * self._number_density(T) * (hbar/self._sTau)


    def _source_electron_0(self, T):
        return self._sBRee * self._number_density(T) * (hbar/self._sTau)
        
    def _source_electron_1(self, T):
        return 0.

    def _source_photon_c(self, E, T):
        EV = self._sE0
        MV = self._sMphi

        x = E/EV

        # If the vector boson mass is smaller than 1 MeV
        # return 0
        if MV < 1.:
            return 0.

        #FSR from V -> e+ e-
        y = me2/MV**2.
        if 1. - 4.*y < x: #Correction to the threshold
            FSRee = 0.
        else:
            W = sqrt(1.-4.*y/(1.-x))
            FSRee = (1./EV) * (alpha/pi) * ( ( 1. + (1.-x)**2. - 4.*y*(x+2.*y)) * log((1.+W)/(1.-W))\
                - ( 1. + (1.-x)**2. + 4.*y*(1.-x))*W )/(1+2.*y)/sqrt(1-4.*y)/x #check the sign 

        #FSR from V -> mu+ mu-
        y = mmu2/MV**2.
        if 1. - 4.*y < x: #Correction to the threshold
            FSRmumu = 0.
        else:
            W = sqrt(1.-4.*y/(1-x))
            FSRmumu = (1./EV) * (alpha/pi) * ( ( 1. + (1.-x)**2. - 4.*y*(x+2.*y)) * log((1.+W)/(1.-W))\
                - ( 1. + (1.-x)**2. + 4.*y*(1.-x))*W )/(1+2.*y)/sqrt(1-4.*y)/x #check the sign 

        #Radiation from mu decay in V -> mu+ mu-
        if ( MV < 2.*mmu or self._sBRmumu < 1.e-15): #Don't calculate if kinematically forbidden or branching ratio is very small/zero
            spec_gamma_mu = 0.
        else:
            #Calculate the Lorentz factors  
            gamma = MV/(2.*mmu)
            beta = sqrt(1.-1./gamma**2.)
            xmu = 2.*E/mmu
        
            #if xmu > (1. - me2/mmu2)*gamma*(1. + beta): #Check the threshold energy
            if E > self._MVspec[0]/2. :
                spec_gamma_mu = 0. 
            else:
                Elog = log10(E) #convert to log10(E) for interpolation
                spec_gamma_mu = 10.**self._spec_gamma[0](Elog)

        #FSR from V -> pi+ pi-
        y = mpic2/MV**2.
        if 1. - 4.*y < x: #Correction to the threshold
            FSRpipi = 0.
        else:
            W = sqrt(1.-4.*y/(1-x)) # Formula in 2003.02273 is wrong
            FSRpipi = (1./EV) * (alpha/pi) * 2.*( (1.-4.*y)*(1-x-2.*y)*log((1.+W)/(1.-W))\
                - (1.-x-x**2.-4.*y*(1.-x))*W )/(1-4.*y)**(3./2.)/x 
            

        #Intermediate FSR (from pion and muon decays) from V -> pi+ pi- 
        #Don't calculate if kinematically forbidden or branching ratio is very small/zero
        if ( MV < 2.*mpic or self._sBRpipi < 1.e-15 ): 
            spec_gamma_pipi = 0.
        elif E > self._MVspec[1]/2. :
            spec_gamma_pipi = 0.
        else:
            #print(E,MV)
            Elog = log10(E) #convert to log10(E) for interpolation
            spec_gamma_pipi = 10.**self._spec_gamma[1](Elog)
        
        
        #For V -> pi0 gamma
        E0 = MV/2.*(1.+mpi2/MV**2.) #Energy of pion
        #Don't calculate if kinematically forbidden or branching ratio is very small/zero
        if ( E0 < mpi or self._sBRpia < 1.e-15 ): 
            spec_gamma_pia = 0.
        else:
            #Egamma = MV/2.*(1.-mpi2/MV**2.) #Energy of photon
            gamma = E0/mpi
            beta = sqrt(1.-1./gamma**2.)
            if E < (1.-beta)*E0/2. or E > (1.+beta)*E0/2.: #Check the threshold energy
                spec_gamma_pia = 0.
            else:
                spec_gamma_pia = 2./(beta*gamma*mpi) 
        
        
        #Gamma ray from V -> pi0 pi+ pi- 
        #Don't calculate if kinematically forbidden or branching ratio is very small/zero
        #Be careful of the threshold since in the spectrum we have taken mpi = mpic
        if ( MV < 2.*mpic + mpi or self._sBR3pi < 1.e-15 ): 
            spec_gamma_3pi = 0.
        elif E > self._MVspec[2]/2. :
            spec_gamma_3pi = 0.
        else:
            Elog = log10(E) #convert to log10(E) for interpolation
            spec_gamma_3pi = 10.**self._spec_gamma[2](Elog)


        gamma_ee = self._sBRee * FSRee 
        gamma_mumu = self._sBRmumu * ( FSRmumu + spec_gamma_mu ) 
        gamma_pipi = self._sBRpipi * ( FSRpipi + spec_gamma_pipi ) 
        gamma_pia = self._sBRpia * spec_gamma_pia 
        gamma_3pi = self._sBR3pi * spec_gamma_3pi 
 
        _sp = self._number_density(T) * (hbar/self._sTau)

        return _sp * ( gamma_ee + gamma_mumu + gamma_pipi + gamma_pia + gamma_3pi )


    def _source_electron_c(self, E, T):
        MV = self._sMphi

        # If the vector boson mass is smaller than 1 MeV
        # return 0
        if MV < 1.:
            return 0.

        #electron spectrum from V -> mu mu  (same for positron)
        #Don't calculate if branching ratio is very small or zero
        if ( MV < 2.*mmu or self._sBRmumu < 1.e-15 ): 
            spec_e_mumu = 0.
        else:
            #Calculate the Lorentz factors    
            gamma = MV/(2.*mmu)
            beta = sqrt(1.-1./gamma**2.)
            xmu = 2.*E/mmu
   
            #if xmu > gamma*(1. + beta*sqrt(1 - 4.*me2/mmu2)):
            if E > self._MVspec[0]/2. :
                spec_e_mumu = 0. 
            else:
                Elog = log10(E) #convert to log10(E) for interpolation
                spec_e_mumu = 10.**self._spec_e[0](Elog)
     
        #electron spectrum from V -> pi+ pi- 
        #Don't calculate if kinematically forbidden or branching ratio is very small/zero
        if ( MV < 2.*mpic or self._sBRpipi < 1.e-15 ): 
            spec_e_pipi = 0.
        else:
            #Calculate the Lorentz factors    
            gammap = mpic/(2.*mmu)*(1. + mmu2/mpic2)
            betap = sqrt(1.-1./gammap**2.)
            xmupmax = gammap*(1. + betap*sqrt(1 - 4.*me2/mmu2))
            
            gamma = MV/(2.*mpic)
            beta = sqrt(1.-1./gamma**2.)
            xmu = 2.*E/mmu
            
            #if xmu > gamma*xmupmax*(1. + beta):
            if E > self._MVspec[1]/2. :
                spec_e_pipi = 0.
            else:
                Elog = log10(E) #convert to log10(E) for interpolation
                spec_e_pipi = 10.**self._spec_e[1](Elog)
        
        #electron spectrum from from V -> pi0 pi+ pi- 
        #Don't calculate if kinematically forbidden or branching ratio is very small/zero
        #Be careful of the threshold since in the spectrum we have taken mpi = mpic
        if ( MV < 2.*mpic + mpi or self._sBR3pi < 1.e-15 ): 
            spec_e_3pi = 0.
        elif E > self._MVspec[2]/2. :
            spec_e_3pi = 0.
        else:
            Elog = log10(E) #convert to log10(E) for interpolation
            spec_e_3pi = 10.**self._spec_e[2](Elog)
     
        
        e_mumu = self._sBRmumu * spec_e_mumu
        e_pipi = self._sBRpipi * spec_e_pipi
        e_3pi = self._sBR3pi * spec_e_3pi 

        _sp = self._number_density(T) * (hbar/self._sTau)

        return _sp * ( e_mumu + e_pipi + e_3pi )


