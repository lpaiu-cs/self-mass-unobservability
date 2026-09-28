# SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later

#coding:utf8

# distutils: language = c++
#-I/usr/share/tempo2/include/ -L/usr/share/tempo2/lib/  -ltempo2 "]
# COmpile avec : make libFittriplecpp.so && python setup_cppFittriple.py build_ext --inplace
# Attention à vérifier que make library fait bien une librairie statique !

#/* This code implements a timing model and fitting procedures aimed at studiyng th triple system J0337+1715 (Ransom et al. 2013)
 #*
 #* Written by Guillaume Voisin 2017 , LUTh, Observatoire de Paris, PSL Research University (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)
 #*/


import numpy as np
import ctypes
from libcpp.vector cimport vector
from libc.stdlib cimport malloc, free

cdef extern from "Diagnostics.h" :
    void IntegralePrems3_1PN_double ( const double s0[6], const double s1[6], const double s2[6],
                         double m0, double m1, double m2,
                         double  center_of_mass[3], double  impulsion[3],
                         double &  energy  )

cdef extern from "Fittriple.h"  : #namespace "shapes":
    cdef cppclass Fittriple:
        Fittriple(char * , char * ) except +
        double Compute_lnposterior(int )

        void Initialise_parameters()

        void Set_fitted_parameter_relativeshifts(vector[double] relativeshift)

        #void Set_parameters_relativeshift(double , double ,
                                      #double , double , double , double , double , double ,
                                      #double , double , double , double , double , double ,
                                      #double , double , double , double )
        #void Set_parameters_relativeshift(double , double ,
                                      #double , double , double , double , double , double ,
                                      #double , double , double , double , double , double ,
                                      #double , double , double , double, double )
        #void Set_parameters_relativeshift(double , double ,
                                      #double , double , double , double , double , double ,
                                      #double , double , double , double , double , double ,
                                      #double , double , double , double, double,
                                      #double)
        #void Set_parameters_relativeshift(double , double ,
                                      #double , double , double , double , double , double ,
                                      #double , double , double , double , double , double ,
                                      #double , double , double , double, double,
                                      #double, double, double )

        void Compute_turns_from_parameters()
        void Create_mask(double threshold)

        void Save_parfile(char *)
        void Save_data(char * filename)
        void Save_output_timing_data(char * filename)
        void Save_masked_copy_of_tim_file(char * masked_tim_file, int * mask)
        void Print()

        void Set_tracker_on()
        void Set_tracker_off()
        int Get_tracker_status()

        void Estimate_parameter_shift_scales(double chi2variation)


        void Set_aberration_delay_on()
        void Set_aberration_delay_off()
        void Set_geometric_delay_on()
        void Set_geometric_delay_off()
        void Set_einstein_delay_on()
        void Set_einstein_delay_off()
        void Set_shapiro_delay_on()
        void Set_shapiro_delay_off()

        void Set_theory(int theorynumber)
        int Get_theory()
        
        void Get_t2parfile(char * t2parfile)
        void Get_specialcase(char * specialcase)

        long int Get_number_of_toas()
        double Get_toa_double(long int toa_number)
        double Get_timeshift_double()
        double Get_treference_double()
        double Get_posepoch_double()
        double Get_residual_double(long int toa_number)
        double Get_error_double(long int toa_number)
        double Get_sat_double(long int toa_number)
        double Get_frequency(long int toa_number)
        double Get_delays(int delaynb, long int toa_number)
        int Get_parameter_set()


        #double* Get_interp_delay(int delaynb)

        long int Get_turn(long int toa_number)

        int Get_n_absparameters()
        int Get_absolute_parameter_map(int param_number)

        int Get_nfitparams()
        double Get_parameter_double(int param_number )
        double Get_reference_parameter_double(int param_number )
        double Get_parameter_scale_double(int param_number )
        int Get_nfitted_params() #Return the number of fitted parameters
        int Get_list_of_fitted_parameters(int param_number) # Return the indexes of fitted parameters. paramènumber < Get_nfitted_params()
        void Get_parameter_name (int paramnumber, char * name) # return the name of each parameter paramnumber
        int Get_nextra()
        double Get_RA_double()
        double Get_DEC_double()
        
        
        double Get_DMX_range_double(int DMXrange_number)
        int Get_number_of_DMX()

        long int Get_number_of_tinterp()
        double Get_tinterp_double(long int tinterp_number )
        double Get_energy_double(long int tinterp_number)
        double Get_center_of_mass_position_double(long int tinterp_number, int component)
        double Get_center_of_mass_impulsion_double(long int tinterp_number, int component)

        long int Get_number_of_removed_toas()
        long int Get_removed_toa(long int toa_number)

        double Get_mass_double(int bodynumber)
        void Get_spinfreq_double(double & truespinfreq, double & truespinfreq1)

        double Get_pulsar_statevector_double(long int tinterp_nb, int component)
        double Get_inner_statevector_double(long int tinterp_nb, int component)
        double Get_outer_statevector_double(long int tinterp_nb, int component)
        double Get_extra_statevector_double(int extra, long int tinterp_nb, int component)
        void Get_analytical_pulsar_statevector_double(long int tinterp_nb, double * statevector)#int component) # Return the analytical, newtonian, state vector computed from the orbital elements at time tinterp[tinterp_nb] with orbel2statevect
        void Get_analytical_inner_statevector_double(long int tinterp_nb, double * statevector)
        void Compute_initial_state_vectors_double(double ** svs)
        
        void Get_fake_bats_and_delays_interp(long int nfakeBats, long double * fake_BATs, long double * fake_delay_geom, long double * fake_delay_ein,
                                                                                    long double * fake_delay_shap, long double * fake_delay_aber)


       # Diagnostics
        void Compute_integrals_of_motion()
        void Get_global_interpolation_values(long int ntimes, long double * times, double * einstein_interp, double * shapiro_aberration_interp)



def Rotate_SSB_to_PSB(alpha, delta, vi):
    sa = np.sin(alpha)
    ca = np.cos(alpha)
    sd = np.sin(delta)
    cd = np.cos(delta)
    
    vx = vi[0] * (-sa)          +  vi[1] * ca     
    vy = vi[0] * (-ca*sd)       +  vi[1] * (  -sa * sd ) + vi[2] *  cd    
    vz = vi[0] * ca*cd          +  vi[1] * sa*cd         + vi[2] *  sd       

    return np.array([vx,vy,vz])




cdef class PyFittriple:

    cdef Fittriple *thisptr      # hold a C++ instance which we're wrapping

    def __cinit__(self, parfile, datafile):
        '''
            fit = PyFittriple(parfile, datafile)
        '''
        print ("datafile : ", datafile)
        self.thisptr = new Fittriple(bytes(parfile, 'utf-8'), bytes(datafile,'utf-8'))#(ctypes.c_char_p(parfile_c), ctypes.c_char_p(datafile_c)) #(parfile, datafile)#

    #def __cinit__(self):
        #parfile = "parfilemcmc"
        #datafile = "timfilemcmc"
        #print ""
        #print " * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
        #print "********** Using default parfile and timfile ***********"
        #print " * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
        #print ""
        #self.thisptr = new Fittriple(parfile, datafile)

    def __dealloc__(self):
        del self.thisptr

    def Compute_lnposterior(self, int fractional = 0):
        return self.thisptr.Compute_lnposterior(fractional)

    def Initialise_parameters(self):
        self.thisptr.Initialise_parameters()
        return

    def Compute_turns_from_parameters(self):
        self.thisptr.Compute_turns_from_parameters()

    def Set_fitted_parameter_relativeshifts(self, relativeshift):
        self.thisptr.Set_fitted_parameter_relativeshifts(relativeshift)
        return

    def Get_number_of_fitted_parameters(self):
        return self.thisptr.Get_nfitted_params()


    def Save_parfile(self, filename) :
        self.thisptr.Save_parfile(bytes(filename, 'utf-8') )

    def Save_data(self, filename) :
        self.thisptr.Save_data(bytes(filename, 'utf-8'))

    def Save_output_timing_data(self, filename):
        self.thisptr.Save_output_timing_data(bytes(filename, 'utf-8'))
    
    def Save_masked_copy_of_tim_file(self, filename, mask):
        '''
        Save_masked_copy_of_tim_file(self, filename, mask)
        
        Save to "filename" a tim file base on the current sorted tim file (i.e. usually the "timfile-sorted" created at startup)
        a file where all the toas with mask = 1 are commented and left untouched if mask = 0.
        '''
        cdef int * cmask;
        cmask = <int *> malloc(mask.size * sizeof(int));
        for i in range(mask.size):
            cmask[i] = mask[i]
            
        self.thisptr.Save_masked_copy_of_tim_file(bytes(filename, 'utf-8'), cmask)
        
        free(cmask)
        
        return 

    def Print(self):
        self.thisptr.Print()

    def Set_tracker_on(self):
        self.thisptr.Set_tracker_on()

    def Set_tracker_off(self):
        self.thisptr.Set_tracker_off()

    def Estimate_parameter_shift_scales(self, chi2variation = 0.01):
        self.thisptr.Estimate_parameter_shift_scales(chi2variation)

    def Create_mask(self, threshold):
        '''
        Create_mask(threshold)
        
        Internally reset the data by removing all the toas which current residual is larger than threshold.
        '''
        self.thisptr.Create_mask(threshold)

    def Get_toas(self):
        ntoas = self.thisptr.Get_number_of_toas()
        toas = np.ndarray(ntoas, dtype=np.float64)
        for i in range(ntoas) :
            toas[i] = self.thisptr.Get_toa_double(i)
        return toas

    def Get_timeshift(self):
        timeshift = self.thisptr.Get_timeshift_double()
        return timeshift


    def Get_treference(self):
        return self.thisptr.Get_treference_double()

    def Get_posepoch(self):
        return self.thisptr.Get_posepoch_double()


    def Get_time_residuals(self):
        '''
            Return the timing residuals in microseconds.
        '''
        ntoas = self.thisptr.Get_number_of_toas()
        residuals = np.ndarray(ntoas, dtype=np.float64)
        for i in range(ntoas) :
            residuals[i] = self.thisptr.Get_residual_double(i)
        return residuals


    def Get_errors(self):
        '''
          Return the array of uncertainties extracted from tempo2's tim file, in microsec.
        '''
        ntoas = self.thisptr.Get_number_of_toas()
        errors = np.ndarray(ntoas, dtype=np.float64)
        for i in range(ntoas) :
            errors[i] = self.thisptr.Get_error_double(i)
        return errors

    def Get_frequencies(self):
        '''
        Return the array of observation frequency extracted from tempo2's tim file, in MHz.
        '''
        ntoas = self.thisptr.Get_number_of_toas()
        freqs = np.ndarray(ntoas, dtype=np.float64)
        for i in range(ntoas) :
            freqs[i] = self.thisptr.Get_frequency(i)
        return freqs

    def Get_SATs(self):
        '''
          Return the array of Site Arrival Times extracted from tempo2's tim file, in MJD.
          WARNING : These are truncated to double precision !
        '''
        print("Don't forget : in python arrival times are truncated to double precision !")
        ntoas = self.thisptr.Get_number_of_toas()
        sats = np.ndarray(ntoas, dtype=np.float64)
        for i in range(ntoas) :
            sats[i] = self.thisptr.Get_sat_double(i)
        return sats


    def Get_delays(self, delaynb):
        '''
        Return delay number delaynb in seconds for each TOA. 
        (Currently disabled) delay < 0 : return some tempo delays
        delaynb = 0 : return the total delay
        delaynb = 1 : return the geometrical delay
        delaynb = 2 : return the sum of the enabled non-geometrical delays besides einstein (shapiro + aberration)
        delaynb = 3 : return the einstein delay
        -- Decomposition of the geometrical delay --
        delaynb = 10 : roemer delay
        delaynb = 11 : schklovski delay
        delaynb = 12 : scklovski correction delay
        delaynb = 13 : kopeikin orbital delay
        delaynb = 14 : kopeikin annual-orbital delay
        delaynb = 15 : kopeikin parallax-orbital delay
        '''
        ntoas = self.thisptr.Get_number_of_toas()
        delays = np.ndarray(ntoas, dtype=np.float64)
        for i in range(ntoas) :
            delays[i] = self.thisptr.Get_delays(delaynb, i)
        return delays

    #def Get_interp_delay(self, delaynb):
        #ninterp = self.thisptr.Get_number_of_tinterp()
        #delays = np.ndarray(ninterp, dtype=np.float64)
        #cdef double * x = self.thisptr.Get_interp_delay(delaynb)
        #for i in range(ninterp) :
            #delays[i] = x[i] * 86400.
        #return delays

    def Get_turns(self):
        ntoas = self.thisptr.Get_number_of_toas()
        turns = np.ndarray(ntoas, dtype=np.int)
        for i in range(ntoas) :
            turns[i] = self.thisptr.Get_turn( i)
        return turns

    def Get_parameters(self):
        npar = self.thisptr.Get_nfitparams()
        pars = np.ndarray(npar, dtype=np.float64)
        for i in range(npar) :
            pars[i] = self.thisptr.Get_parameter_double(i)
        return pars

    def Get_reference_parameters(self):
        npar = self.thisptr.Get_nfitparams()
        pars = np.ndarray(npar, dtype=np.float64)
        for i in range(npar) :
            pars[i] = self.thisptr.Get_reference_parameter_double(i)
        return pars

    def Get_parameter_scales( self ) :
        npar = self.thisptr.Get_nfitparams()
        parscale = np.ndarray(npar, dtype=np.float64)
        for i in range(npar) :
            parscale[i] = self.thisptr.Get_parameter_scale_double(i)
        return parscale

    def Get_fitted_parameter_map(self) :
        '''
            Return the map of the indexes of the fitted parameters . Same as Get_fitted_parameters_list, but not sorted;
        '''
        nparfitted = self.thisptr.Get_nfitted_params() #Return the number of fitted parameters
        fittedlist = np.ndarray(nparfitted, dtype=np.int)
        for i in range(nparfitted):
            fittedlist[i] = self.thisptr.Get_list_of_fitted_parameters(i)
        return fittedlist

    def Get_absolute_parameter_map(self):
        nabsparam = self.thisptr.Get_n_absparameters()
        absparams = np.ndarray(nabsparam, dtype=np.int)
        for i in range(nabsparam):
            absparams[i] = self.thisptr.Get_absolute_parameter_map(i)
        return absparams
        #cdef int * cabsparams
        #cabsparams = self.thisptr.Get_absolute_parameter_map()
        #for i in range(nabsparam):
            #absparams[i] = cabsparams[i]
        #free(cabsparams)
        #return absparams

    def Get_fitted_parameters_list(self) :
        '''
            Return the sorted list of the indexes of the fitted parameters in parameters = Get_parameters()
        '''
        nparfitted = self.thisptr.Get_nfitted_params() #Return the number of fitted parameters
        fittedlist = np.ndarray(nparfitted, dtype=np.int)
        for i in range(nparfitted):
            fittedlist[i] = self.thisptr.Get_list_of_fitted_parameters(i)
        fittedlist.sort()
        return fittedlist

    def Get_parameter_names (self) :
        '''
            Return the list of the names of all the parameters.
        '''
        cdef char onename[100]
        npar = self.thisptr.Get_nfitparams()
        parnames = []
        print(npar)
        for i in range(npar) :
             self.thisptr.Get_parameter_name(i, onename)
             parnames.append(str(onename, 'UTF-8')) 
        return parnames

    def Get_DMX_ranges(self):
        ndmx = self.thisptr.Get_number_of_DMX()
        DMXr = np.zeros(ndmx+1, dtype=np.float64)
        for i in range(DMXr.size):
            DMXr[i] = self.thisptr.Get_DMX_range_double(i)
        return DMXr

    def Compute_integrals_of_motion(self) :
        '''
            Need to be run BEFORE getting any integral of motion.
            Need to be run AFTER Compute_lnposterior()
        '''
        self.thisptr.Compute_integrals_of_motion()

    def Get_energies(self) :
        ninterp = self.thisptr.Get_number_of_tinterp()
        x = np.ndarray(ninterp, dtype=np.float64)
        for i in range(ninterp) :
            x[i] = self.thisptr.Get_energy_double(i)
        return x

    def Get_center_of_mass_positions(self) :
        ninterp = self.thisptr.Get_number_of_tinterp()
        x = np.ndarray((ninterp,3), dtype=np.float64)
        for i in range(ninterp) :
            for j in range(3) :
                x[i,j] = self.thisptr.Get_center_of_mass_position_double(i,j)
        return x

    def Get_center_of_mass_impulsions(self) :
        ninterp = self.thisptr.Get_number_of_tinterp()
        x = np.ndarray((ninterp,3), dtype=np.float64)
        for i in range(ninterp) :
            for j in range(3) :
                x[i,j] = self.thisptr.Get_center_of_mass_impulsion_double(i,j)
        return x

    def Get_tinterp(self) :
        ninterp = self.thisptr.Get_number_of_tinterp()
        x = np.ndarray(ninterp, dtype=np.float64)
        for i in range(ninterp) :
            x[i] = self.thisptr.Get_tinterp_double(i)
        return x


    def Get_removed_toas(self) :
        n = self.thisptr.Get_number_of_removed_toas()
        x = np.ndarray(n, dtype=np.int)
        for i in range(n) :
            x[i] = self.thisptr.Get_removed_toa(i)
        return x

    def Get_nextra(self):
        '''
            Return number of extra bodies (beyond the 3 basic ones)
        '''
        return self.thisptr.Get_nextra()
    
    def Get_masses(self) :
        '''
            Return masses
        '''
        nbody = self.Get_nextra()+3
        Ms = np.zeros(nbody)
        for i in range(nbody):
            Ms[i] = self.thisptr.Get_mass_double(i)
        return Ms
    #np.array([self.thisptr.Get_mass_double(0), self.thisptr.Get_mass_double(1), self.thisptr.Get_mass_double(2)], dtype=np.float64)

    def Get_true_spinfreq(self):
        cdef double truespinfreq
        cdef double truespinfreq1
        self.thisptr.Get_spinfreq_double(truespinfreq, truespinfreq1)
        return truespinfreq, truespinfreq1

    def Get_interp_state_vectors(self, PSB=False) :
        '''
        state_vector_pulsar, state_vector_inner_companion, state_vector_outer_companion [, state_vector_extra_companions]  = Get_interp_state_vectors()

        Return the state vectors for each body at every interpolation times. If no extra body is present then return only 3 arrays, otherwise return an extra array with shape = (nextra, ninterp, 6)
        
        PSB : if False then return vectors in the SSB reference frame, otherwise rotate them such that the z coordinate is along the line of sight
        '''
        ninterp = self.thisptr.Get_number_of_tinterp()
        sp = np.ndarray((ninterp,6), dtype=np.float64)
        si = np.ndarray((ninterp,6), dtype=np.float64)
        so = np.ndarray((ninterp,6), dtype=np.float64)
        for i in range(ninterp) :
            for j in range(6) :
                sp[i,j] = self.thisptr.Get_pulsar_statevector_double(i, j)
                si[i,j] = self.thisptr.Get_inner_statevector_double(i, j)
                so[i,j] = self.thisptr.Get_outer_statevector_double(i, j)
        
        if PSB :
            RA = self.thisptr.Get_RA_double()
            DEC = self.thisptr.Get_DEC_double()
            for i in range(ninterp) :
                sp[i,:3] = Rotate_SSB_to_PSB(RA,DEC,sp[i,:3])
                si[i,:3] = Rotate_SSB_to_PSB(RA,DEC,si[i,:3])
                so[i,:3] = Rotate_SSB_to_PSB(RA,DEC,so[i,:3])
                sp[i,3:] = Rotate_SSB_to_PSB(RA,DEC,sp[i,3:])
                si[i,3:] = Rotate_SSB_to_PSB(RA,DEC,si[i,3:])
                so[i,3:] = Rotate_SSB_to_PSB(RA,DEC,so[i,3:])
        
        nextra = self.thisptr.Get_nextra()
        if nextra > 0:
            sextra = np.ndarray((nextra, ninterp, 6), dtype=np.float64)
            for k in range(nextra):
                for i in range(ninterp) :
                    for j in range(6) :
                        sextra[k, i,j] = self.thisptr.Get_extra_statevector_double(k, i, j)
                if PSB : 
                    for i in range(ninterp) :
                        sextra[k, i, :3] = Rotate_SSB_to_PSB(RA,DEC,sextra[k, i, :3])
                        sextra[k, i, 3:] = Rotate_SSB_to_PSB(RA,DEC,sextra[k, i, 3:])
                        
            return sp, si, so, sextra
        else : 
            return sp, si, so


    def Get_fake_bats_and_delays_interp(self,n ):
        '''
            bats, delay_geom, delay_ein, delay_shap, delay_aber = fit.Get_fake_bats_and_delays_interp(n)

            Generate "n" fake delays and bats close to each interpolation point except on the edges of the grid.
            The fake toas are computed according to the current parameters.
        '''
        cdef long int nbats = n #self.thisptr.Get_number_of_tinterp() - 2
        cdef long double * cbats =  <long double*>malloc( sizeof(long  double) * nbats )
        cdef long double * cdelay_geom =  <long double*>malloc( sizeof(long  double) * nbats )
        cdef long double * cdelay_ein=  <long double*>malloc( sizeof(long  double) * nbats )
        cdef long double * cdelay_shap =  <long double*>malloc( sizeof(long  double) * nbats )
        cdef long double * cdelay_aber =  <long double*>malloc( sizeof(long  double) * nbats )



        self.thisptr.Get_fake_bats_and_delays_interp(nbats, cbats, cdelay_geom, cdelay_ein, cdelay_shap, cdelay_aber)

        bats = np.ndarray(nbats)
        delay_geom = np.ndarray(nbats)
        delay_ein = np.ndarray(nbats)
        delay_shap = np.ndarray(nbats)
        delay_aber = np.ndarray(nbats)
        for i in range(nbats) :
            delay_geom[i] = cdelay_geom[i]
            delay_ein[i] = cdelay_ein[i]
            delay_shap[i] = cdelay_shap[i]
            delay_aber[i] = cdelay_aber[i]
            bats[i] = cbats[i]

        free(cdelay_geom)
        free(cdelay_ein)
        free(cdelay_shap)
        free(cdelay_aber)
        free(cbats)

        return bats, delay_geom, delay_ein, delay_shap, delay_aber


    def Get_global_interpolation_values(self, times):
        cdef long int ntimes = times.size
        cdef long double * ctimes = <long double*>malloc( sizeof(long  double) * ntimes )
        cdef double * einstein_interp = <double*>malloc( sizeof(double) * ntimes )
        cdef double * shapiro_aberration_interp = <double*>malloc( sizeof(double) * ntimes )

        for i in range(ntimes) :
            ctimes[i] = times[i]

        self.thisptr.Get_global_interpolation_values(ntimes, ctimes, einstein_interp, shapiro_aberration_interp)

        np_einstein_interp = np.ndarray(ntimes)
        np_shapiro_aberration_interp = np.ndarray(ntimes)
        for i in range(ntimes) :
            np_einstein_interp[i] = einstein_interp[i]
            np_shapiro_aberration_interp[i] = shapiro_aberration_interp[i]

        free(einstein_interp)
        free(shapiro_aberration_interp)

        return np_einstein_interp, np_shapiro_aberration_interp


    def Get_analytical_interp_state_vectors(self) :
        '''
             Return the analytical, newtonian, state vector computed from the orbital elements at time tinterp[tinterp_nb] with orbel2statevect
        '''
        ninterp = self.thisptr.Get_number_of_tinterp()
        sp = np.ndarray((ninterp,6), dtype=np.float64)
        cdef double spa[6]
        cdef double sia[6]
        cdef double szero[6]
        #cdef double center_of_mass[3]
    #    cdef double impulsion[3]
        #cdef double energy
        si = np.ndarray((ninterp,6), dtype=np.float64)
        #cdm = np.ndarray((ninterp,3), dtype=np.float64)
    #    imp = np.ndarray((ninterp,3), dtype=np.float64)
    #    nrj = np.ndarray(ninterp, dtype=np.float64)
        #so = np.ndarray((ninterp,6), dtype=np.float64)
        for i in range(3) :
            szero[i] = 0.

        Mp, Mi, Mo = self.Get_masses()

        for i in range(ninterp) :
            self.thisptr.Get_analytical_pulsar_statevector_double(i, spa)
            self.thisptr.Get_analytical_inner_statevector_double(i, sia)
            # IntegralePrems3_1PN_double ( spa, sia, szero,
            #              Mp, Mi,  Mo,
            #              center_of_mass, impulsion,
            #              energy  )
            for j in range(6) :
                sp[i,j] = spa[j]
                si[i,j] = sia[j]
            #for j in range(3) :
            #    cdm[i,j] = center_of_mass[j]
            #    imp[i,j] = impulsion[j]
            #nrj[i] = energy
        #    #for j in range(6) :
        #        #sp[i,j] = self.thisptr.Get_analytical_pulsar_statevector_double(i, j)
        #        ##si[i,j] = self.thisptr.Get_pulsar_statevector_double(i, j)
        #        #sp[i,j] = self.thisptr.Get_pulsar_statevector_double(i, j)
        return sp, si#, cdm, imp, nrj

    def Compute_initial_state_vectors(self):
        '''
            Compute the initial state vectors according to current parameters 
            and return them under the form state_vectors[nbody,6] where the 
            first 3 components are position and the last 3 are velocity. 
            
            This routine allows to recompute only the initial state vectors without having to integrate over the whole time span as well (using self.Compute_lnposterior()).
            
            It also recomputes all the corresponding derived quantities, in particular masses.
        '''
        cdef int nbody = 3+self.thisptr.Get_nextra()
        statevects = np.zeros((nbody, 6))
        cdef double ** svs
        svs = <double **> malloc(sizeof(double*) * nbody)
        for i in range(nbody):
            svs[i] = <double *> malloc(sizeof(double) * 6)
            
        self.thisptr.Compute_initial_state_vectors_double(svs)
        
        for i in range(nbody):
            for j in range(6):
                statevects[i,j] = svs[i][j]
        
        for i in range(nbody):
            free(svs[i])
        free(svs)
        
        
        return statevects
    
    
    def  Set_aberration_delay_on(self):
        self.thisptr.Set_aberration_delay_on()
        return

    def  Set_aberration_delay_off(self):
        self.thisptr.Set_aberration_delay_off()
        return

    def  Set_geometric_delay_on(self):
        self.thisptr.Set_geometric_delay_on()
        return

    def  Set_geometric_delay_off(self):
        self.thisptr.Set_geometric_delay_off()
        return

    def  Set_einstein_delay_on(self):
        self.thisptr.Set_einstein_delay_on()
        return

    def  Set_einstein_delay_off(self):
        self.thisptr.Set_einstein_delay_off()
        return

    def  Set_shapiro_delay_on(self):
        self.thisptr.Set_shapiro_delay_on()
        return

    def  Set_shapiro_delay_off(self):
        self.thisptr.Set_shapiro_delay_off()
        return

    def Set_theory(self, theorynumber):
        self.thisptr.Set_theory(theorynumber)
        return

    def Get_theory(self):
        return self.thisptr.Get_theory()

    def Get_tracker_status(self):
        return self.thisptr.Get_tracker_status()

    def Get_parameter_set(self):
        return self.thisptr.Get_parameter_set()
    
    def Get_t2parfile(self):
        cdef char t2parfile[500]
        return str(self.thisptr.Get_t2parfile(t2parfile), 'UTF-8')
    
    def Get_specialcase(self):
        cdef char specialcase[500]
        return str(self.thisptr.Get_specialcase(specialcase), 'UTF-8')
    
