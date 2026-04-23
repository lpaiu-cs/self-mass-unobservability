# SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later

# This file contains tools to plot the results of integrate.exe

import numpy as np
from os.path import exists
import matplotlib.pyplot as pplot
import celestial_mechanics_tools as cm
import Constantes_physiques as cp

class integration:
    def __init__(self, fileprefix="trajectories",Ms=None, dir_obs='reference'):     
        '''
        fileprefix : prefix of filenames to load. 
                - Trajectory file : "fileprefix"
                - Integral of motion : "fileprefix"+"-intofmotion"
        
        Ms : numpy array containing the masses of each bodies in Msun. Necessary for computing orbital elements
        
        dir_obs : defines the plane with respect to which space orientation is defined for computation of orbital elements:
                - 'reference' : take the reference plane with normal = total angular momentum
                - 'timing' : normal is the z axis 
                - array-like : vector of the normal direction
        '''
        trajfl = np.loadtxt(fileprefix)
        self.ts = trajfl[:,0]
        self.nt = self.ts.size
        self.nbody = np.int((trajfl.shape[1] -1 )/6)
        self.xs = np.zeros((self.nbody, self.nt, 3))
        self.vs = np.zeros((self.nbody, self.nt, 3))
        for i in range(self.nbody):
            self.xs[i,:,:] = trajfl[:,(i*6 +1):(i*6+3+1)]
            self.vs[i,:,:] = trajfl[:,(i*6 +3+1):(i*6+6+1)]
            
        if exists(fileprefix+"-intofmotion"):
            imofl = np.loadtxt(fileprefix+"-intofmotion")
            self.flag_imo =True
            self.ecof = imofl[:,1]
            self.xcof = imofl[:, 2:5]
            self.vcof = imofl[:, 5:8]
            self.L = imofl[:, 8:11]
        else : 
            self.flag_imofl=False
            
        if Ms is not None:
            self.Ms = Ms
        else:
            self.Ms = np.zeros(self.nbody)
            print("Warning : you need to give masses to be able to compute orbital elements !")
            
        self.dir_obs = dir_obs
            
        return 
    
    def Plot_integrals_of_motion(self, figsize =(10,8)):
        if self.flag_imo==False:
            return
        
        figx = pplot.figure(figsize=figsize)
        plt = figx.add_subplot(111)
        plt.plot(self.ts, self.xcof[:,0], label='$x$')
        plt.plot(self.ts, self.xcof[:,1], label='$y$')
        plt.plot(self.ts, self.xcof[:,2], label='$z$')
        plt.set_xlabel("t (days)")
        plt.set_ylabel("Position (m)")
        pplot.legend()
        figx.show()
        
        figv = pplot.figure(figsize=figsize)
        plt = figv.add_subplot(111)
        plt.plot(self.ts, self.vcof[:,0], label='$v_x$')
        plt.plot(self.ts, self.vcof[:,1], label='$v_y$')
        plt.plot(self.ts, self.vcof[:,2], label='$v_z$')
        plt.set_xlabel("t (days)")
        plt.set_ylabel("Velocity (m/s)")
        pplot.legend()
        figv.show()
        
        fige = pplot.figure(figsize=figsize)
        plt = fige.add_subplot(111)
        plt.plot(self.ts, (self.ecof - self.ecof[0])/self.ecof[0], label=r'$\Delta E/E_0$')
        plt.set_xlabel("t (days)")
        plt.set_ylabel("Relative variation")
        pplot.legend()
        fige.show()
        
        figL = pplot.figure(figsize=figsize)
        plt = figL.add_subplot(111)
        plt.plot(self.ts, (self.L[:,0] - self.L[0,0])/self.L[0,0], label=r'$\Delta L_x/L_x^0$')
        plt.plot(self.ts, (self.L[:,1] - self.L[0,1])/self.L[0,1], label=r'$\Delta L_y/L_y^0$')
        plt.plot(self.ts, (self.L[:,2] - self.L[0,2])/self.L[0,2], label=r'$\Delta L_z/L_z^0$')
        plt.set_xlabel("t (days)")
        plt.set_ylabel("Relative variation")
        pplot.legend()
        figL.show()
        
        return figx, figv, fige, figL
    
    def Plot_components(self, figsize =(10,8)):
        
        for i in range(self.nbody):
            fig = pplot.figure(figsize=figsize)
            plt = fig.add_subplot(111)
            plt.plot(self.ts, self.xs[i,:,0], label='$x_{:}$(m)'.format(i))
            plt.plot(self.ts, self.xs[i,:,1], label='$y_{:}$(m)'.format(i))
            plt.plot(self.ts, self.xs[i,:,2], label='$z_{:}$(m)'.format(i))
            plt.set_xlabel("t (days)")
            plt.set_ylabel("Position (m)")
            pplot.legend()
            pltv= plt.twinx()
            pltv.plot(self.ts, self.vs[i,:,0], label='$vx_{:}$(m/s)'.format(i), linestyle='--')
            pltv.plot(self.ts, self.vs[i,:,1], label='$vy_{:}$(m/s)'.format(i), linestyle='--')
            pltv.plot(self.ts, self.vs[i,:,2], label='$vz_{:}$(m/s)'.format(i), linestyle='--')
            pltv.set_ylabel("Velocity (m/s)")
            pplot.legend()
            fig.show()
                
        return 
    
    
    def Compute_orbital_elements(self):
        '''
        Compute a list of dictionaries of orbital elements and other derived quantities at each times and for each subset of bodies in the hierarchical system. 
        For example, return orbel such that:
            - orbel[0]['aorb'][3] = separation between body 1 and body 0 at time self.ts[3]. 
            - orbel[2]['Porb'][4] = orbital period corresponding to the orbit of body 3 (i.e. 4th body) with the centre of mass of bodies 0, 1 and 2 at time self.ts[4].
            
        
        '''
        #if Ms is None:
            #if self.Ms is None : 
                #print("You need to provide Ms if you want orbital elements to be computed!")
                #return
            #else:
                #Ms = self.Ms
                #print("ms ssss ", Ms)
         
        Ms = self.Ms
        
        if self.dir_obs == 'reference' :
            dirobs = 0
            for i in range(self.nbody):
                dirobs += np.cross(self.xs[i,0,:], self.vs[i,0,:]) * Ms[i]
        elif self.dir_obs == "timing":
            dirobs = np.array([0,0,1])
        elif type(self.dir_obs) is not str :
            dirobs = np.array(self.dir_obs)
        
        dirobs /= np.linalg.norm(dirobs)
            
        print("Plane = ", self.dir_obs, "Normal to plane = ", dirobs)
        
        s2 = np.zeros(6)
        s1 = np.zeros(6)
        
        # Get the keys
        s1 = xv2sv(self.xs[0, 0,:], self.vs[0, 0,:])
        mu=cp.Ggrav * Ms[:1].sum()
        orbel = State_vectors_to_orbital_elements(s1, mu, dir_obs=dirobs, t=self.ts[0])
        keys = list(orbel.keys())
        shapes=[]
        for k in keys:
            shape = np.zeros(len(orbel[k].shape) +1, dtype=int)
            shape[0] = self.nt
            shape[1:] = orbel[k].shape
            shapes.append(shape)
            
        self.orbels = []
        # Fill dictionary
        for nb in range(1, self.nbody):
            mu = cp.Ggrav * Ms[:nb+1].sum()
            orbel = {}
            for i in range(len(keys)):
                orbel[keys[i]] = np.zeros(shapes[i])
            for i in range(self.nt):
                s1 = xv2sv(np.average(self.xs[:nb,i,:], axis = 0, weights=Ms[:nb]), np.average(self.vs[:nb,i,:], axis = 0, weights= Ms[:nb]))
                s2 = xv2sv(self.xs[nb, i,:], self.vs[nb, i,:])
                dico = State_vectors_to_orbital_elements(s2-s1, mu, dir_obs=dirobs, t=self.ts[i])
                for k in keys:
                    orbel[k][i] = dico[k]
            self.orbels.append(orbel)
                
        return self.orbels.copy()
    
    def __check_orbels(self):
        '''
        Check if orbital elements have already been computed and if not compute them.
        Return the corresponding dictionary
        '''
        try:
            orbels = self.orbels
            print("Using already computed orbital elements")
        except :
            print("Computing orbital elements...")
            orbels = self.Compute_orbital_elements()
        return orbels
    
    
    
    def Get_orbels(self, time_index):
        orbels = self.__check_orbels()
        orbels = []
        print("*** Masses (Msun):")
        for i in range(self.nbody):
            print("{:.10}".format(self.Ms[i]/cp.Msol))
        print()
        for i in range(self.nbody-1):
            print("*** Subset {:} :".format( i ) )
            dico = {}
            for k in self.orbels[i].keys():
                print(k + "  ", self.orbels[i][k][time_index])
                dico[k] = self.orbels[i][k][time_index]
            orbels.append(dico)
            print()
        return orbels
    
    
    
    def Plot_orbital_elements(self, subset, elements, figsize = (10,8)):
        '''
            Plot the orbital elements in "elements" of "subset" vs time. 
            
            Parameters:
                * subset : integer from 0 to N-1, where N is the number of bodies, giving the corresponding subset in the hierarchy of systems (see also Compute_orbital_elements)
                * elements : list of strings containing the names orbital elements or derived quantities to plot. Any of : ['aorb', 'Porb', 'norb', 'ecc', 'incl', 'omperi', 'tperi', 'Oman', tasc_approx', 'true_anomaly','mean_anomaly','eccentricity_anomaly', 'mean_anomaly', 'Laplace_vector', 'ang_momentum', 'line_ascnode']
            
        '''
        orbels = self.__check_orbels()
        
        for k in elements:
            fig = pplot.figure(figsize=figsize)
            plt = fig.add_subplot(111)
            plt.set_title('{:}: '.format(subset) + k)
            plt.plot(self.ts, orbels[subset][k])
            fig.show()
        return 
    
    
    def Plot_orbital_element_map(self, subset1, element1, subset2, element2, 
                                 figsize = (10,8), markersize=4, alpha=0.5):
        '''
            Plot the orbital elements in "elements" of "subset" vs time. 
            
            Parameters:
                * subset : integer from 0 to N-1, where N is the number of bodies, giving the corresponding subset in the hierarchy of systems (see also Compute_orbital_elements)
                * elements : list of strings containing the names orbital elements or derived quantities to plot. Any of : ['aorb', 'Porb', 'norb', 'ecc', 'incl', 'omperi', 'tperi', 'Oman', tasc_approx', 'true_anomaly','mean_anomaly','eccentricity_anomaly', 'mean_anomaly', 'Laplace_vector', 'ang_momentum', 'line_ascnode']
                * Ms : array of masses corresponding to the 4 bodies (in principle could be derived from xs and vs but too complicated to be worth). 
            
        '''
        orbels = self.__check_orbels()
        
        fig = pplot.figure(figsize=figsize)
        plt = fig.add_subplot(111)
        plt.scatter(orbels[subset1][element1], orbels[subset2][element2], s=markersize, alpha=alpha)
        plt.set_ylabel('{:}: '.format(subset2) + element2)
        plt.set_xlabel('{:}: '.format(subset1) + element1)
        fig.show()
        return 
        
        
def xv2sv(x,v):
    '''
    return the state vector corresponding to tposition x and velocity v
    '''
    sv = np.zeros(6)
    sv[0:3] = x
    sv[3:6] = v
    return sv
        
def Derive_mu_from_2_sv(sv1,sv2):
    '''
        mu = GMtot can in principle be computed from x and v, but I don't know a simple formula.. 
        On the other hand if one has two different state vectors sv1 and sv2 from different instants 
        (but exactly not a period apart)  then mu is the coefficient such that q1=q2 where q is the Laplace vector. 
        Beware, this assumes there is no perturbations i.e. q is constant
    '''
    r1 = sv1[:3]
    nr1 = np.linalg.norm(r1)
    rn1 = r1 / nr1
    v1 = sv1[3:6] #- s2[3:6]
    h = np.cross(r1,v1) # Angular Momentum
    
    r2 = sv2[:3]
    nr2 = np.linalg.norm(r2)
    rn2 = r2 / nr2
    v2 = sv2[3:6] 
       
    drdot_x_h = np.cross(v1 - v2,h)
    mus = drdot_x_h/(rn1 - rn2)
    
    return mus


def Derive_mu_from_sv(sv, maxit=100, eps=1e-10):
    '''
        mu = GMtot 
    '''
    r = sv[:3]
    nr = np.linalg.norm(r)
    rn = r / nr
    v = sv[3:6]
    nv = np.linalg.norm(v)
    h = np.cross(r,v) # Angular Momentum
    h2 = h[0]**2 + h[1]**2 + h[2]**2
    rdot_x_h = np.cross(v,h)
    
    mu0 = nv**2*nr
    h2_r = h2/nr
    mu = mu0 *0.8
    oldmu=0.
    i=0
    print(mu, oldmu, i, maxit,(abs(mu - oldmu)/mu))
    while ((abs(mu - oldmu)/mu > eps) and (i < maxit)) :
        qlap = rdot_x_h - mu * rn # Laplace vector
        nqlap = np.linalg.norm(qlap)
        ecc = nqlap / mu
        oldmu = mu 
        mu = (mu0 - h2_r /(1.+ecc))/ecc
        i+=1
        print(i, mu, oldmu, ecc, mu0-h2_r/(1+ecc)) 
    
    if i == maxit : 
        print("Warning: maxit reached in Derive_mu_from_sv !")
        
    return mu
    

def State_vectors_to_orbital_elements(sv, mu, dir_obs=[0.,0.,1.], t=0) : 
    '''
        sv : [x,y,z, vx, vy, vz] in meters and meter/s is the relative state vector (ie state vector seen from one of the bodies)
        dir_obs : vector normal to the plane of the sky
        t : time of observation of sv (in seconds)
        mu : GMtot in SI units
    '''
    kobs = np.array(dir_obs) / np.linalg.norm(dir_obs) # vector orthogonal to plane of sky
    r = sv[:3] 
    nr = np.linalg.norm(r)
    rn = r / nr
    v = sv[3:6] 
    h = np.cross(r,v) # Angular Momentum
    nh = np.linalg.norm(h)
    h2 = nh**2
    nv = np.linalg.norm(v)
    rdot2 = nv**2
    rdot_x_h = np.cross(v,h)

    


    qlap = rdot_x_h - mu * r/nr # Laplace vector
    nqlap = np.linalg.norm(qlap)
    easc = np.cross(kobs, h) # unit vector of the line of ascending nodes
    easc /= np.linalg.norm(easc)

    # Defining a frame (xobs, yobs, kobs) by projection of x,y,z on the reference/observer's plane
    xobs = np.array([1,0,0]) - np.array([1,0,0]).dot(kobs)*kobs # projection x base vector on the plane of the observer (or plane of reference) 
    xobs /= np.linalg.norm(xobs)
    yobs = np.array([0,1,0]) - np.array([0,1,0]).dot(kobs)*kobs # projection y base vector on the plane of the observer (or plane of reference) 
    yobs /= np.linalg.norm(yobs)
    
    orientation = np.sign(np.array([0,0,1]).dot(kobs)) # if orientation>0 then orientation of the frame (xobs, yobs, kobs) is the same as that of (x,y,z)
    Oman = orientation*np.arctan2(easc.dot(yobs), easc.dot(xobs)) # Longitude of ascending node

    slr = h2 / mu # semi-latus rectum
    anglei = np.arccos(h.dot(kobs) / nh)
    ecc = nqlap / mu
    sep = slr/ (1. - ecc**2)
    OmPeri = np.arccos(np.dot(easc, qlap)/nqlap)
    Porb = np.pi * 2 * np.sqrt(sep**3/mu)
    MeanMotion = np.pi * 2 / Porb # Mean motion

    # Frame of the orbit (qlap, yorb, h)
    xorb = qlap / nqlap
    yorb = np.cross(h, qlap)
    yorb /= np.linalg.norm(yorb)
    vtrue = np.arctan2(np.dot(r, yorb), np.dot(r, xorb)) # True anomaly
    sinE = np.dot(r, yorb) / (sep * np.sqrt(1 - ecc**2))
    cosE = np.dot(r, xorb)/ sep + ecc
    Eano = np.arctan2(sinE, cosE) # eccentric anomaly
    MeanAnomaly = Eano - sinE * ecc # Mean anomaly
    tperi = t - MeanAnomaly / MeanMotion # Time of passage at periastron
    tasc_approx = tperi - OmPeri / MeanMotion # Time of passage at ascending node to zeroth order in eccentricity

    # Make a dictionary
    dico = {'aorb':sep, 'Porb':Porb, 'norb' : MeanMotion, 'ecc':ecc,
            'incl':anglei, 'omperi':OmPeri, 'tperi': tperi, 'Oman':Oman,
            'tasc_approx':tasc_approx, 'true_anomaly':vtrue,
            'mean_anomaly': MeanAnomaly,
            'eccentricity_anomaly':Eano, 'mean_anomaly':MeanAnomaly,
            'Laplace_vector':qlap, 'ang_momentum':h, 'line_ascnode':easc}

    return dico 


        
