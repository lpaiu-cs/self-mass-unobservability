# SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later

#coding:utf8
import numpy as np
from outils import inversetrigo

# CODATA : P.J. Mohr et B.N. Taylor, Rev. Mod. Phys., 2005, 77,1-107
# !! EN 2018, les unités devraient être redéf par le BIPM

pi = np.arccos(-1.)
deuxpi = 2.*pi
clight = 299792458.	#    ! m/s Par déf
eps0 = 10.**7/(4.*pi*clight**2)  #! F/m
mu0 = 1./(clight**2 * eps0)
hbarre = 1.054571800e-34 #(13) #! constante de Planck/2pi en (Js) (CODATA 2014)
asf = 7.2973525664*10**(-3) #(17)! constante de structure fine (CODATA 2014)
echarge = -1.6021766208e-19 #(98) !/sqrt(eps0*4.*pi); ! C (CODATA 2014) 
emass = 9.10938356e-31 #(11)# kg CODATA 2014
protonmass = 1.672621898e-27 # kg, proton mass (source : CODATA http://physics.nist.gov/cgi-bin/cuu/Value?mp|search_for=proton+mass)
lcompton = hbarre / (emass * clight) # longueur d'onde compton de l'electron


Bcrit = emass**2 * clight**2/(np.abs(echarge)*hbarre)

elcompton = hbarre/(emass*clight) # Longueur d'onde Compton de l'electron
kb = 1.38064852 * 10.**(-23) #(79) Constante de Boltzmann (J/K) (Par déf, CODATA 2014, http://physics.nist.gov/cgi-bin/cuu/Value?k|search_for=boltzman)
Ggrav = 6.67408e-11      # (31) m^3 s^-2 kg^-1 (CODATA 2014)
Msol = 1.988435e30         # Masse du Soleil (kg)
daysec = 86400.				# secondes par jour julien

raieHalpha = 656.3 # nanametre raie Halpha

eVenJoule = abs(echarge) # J / eV
emass_eV = emass * clight**2 / eVenJoule # electron mass in eV


Telectron = emass * clight**2 / kb # energy of the electron in kelvins 

meterly = clight * 86400. * 365.25 # meter per light year
meterua = clight * 8. *60. # meter per astronomical unit
arcsec = 1./3600. * deuxpi / 360. # radian per arcsec
arcsech = deuxpi / (24. * 3600. )
meterpc = meterua / (2.*np.tan(0.5*arcsec))  # meters / parsec


def pulsar_ppdot(p = 1., pdot=10.**(-15), istar = 10.**(39), rstar = 10.**4) :
    Bpole = np.sqrt(pdot *p * istar * clight**3 * mu0) / deuxpi**2 / rstar**3
    Lsd = istar * deuxpi**2 / p**3 *pdot
    rhoc = 2.*eps0*deuxpi/p * Bpole / np.abs(echarge) 
    print('Mag field (dipolar at the pole, in Teslas) %e'%Bpole)
    print('Spin down power (Watts) %e'%Lsd)
    print('Electron/positron density at the surface (m^(-3) %e', rhoc)
    return Bpole, Lsd

def pulsar(distsurf = 10**4, rayon = 10.**4, Bsurf=10.**8, periode=1., incl=0.5*pi ):
    '''
    Donne tout un tas de grandeurs caractéristiques pour un pulsar.
    incl = angle between magnetic and rotation axis
    '''
    B = Bsurf * (rayon/(rayon + distsurf))**3
    freq = 1./periode
    puls = 2.*pi/periode
    rlc = clight / puls
    champE = puls * (distsurf + rayon) * B
    gammamax = (1.5 * 4.* pi * eps0 * champE * (rayon+distsurf)**2 / np.abs(echarge))**(0.25)
    particledensity = 7.*10.**(2) * Bsurf * np.cos(incl) / periode # Goldreich and Julian, 1969 , formula 9
    rlarmor = emass * clight/np.abs(echarge)/Bsurf
    lambdap = clight *np.sqrt(Foremass *eps0 / (echarge**2 * particledensity * 10.**6) ) # Epaisseur de peau = c/ omega_p
    print("Champ mag dipolaire à %e m = %e "%( distsurf, B))
    print("Champ électrique force-free correspondant = %e"%champE)
    print("Facteur de Lorentz max d'un électron freiné par courbure %e"%gammamax)
    print("Fréquence d'un pulsar de période %e s = %e s^-1"%(periode, freq))
    print("Pulsation d'un pulsar de période %e s = %e s^-1"%(periode, puls))
    print("Cylindre de lumière d'un pulsar de période %e s = %e m"%(periode, rlc))
    print("Particle number density near the surface in cm^-3 = %e"%particledensity)
    print("Rayon de larmor à la surface pour une particule avec vperp = c et gammaperp=1 m = %e"%rlarmor)
    print("Epaisseur de peau à la surface in m = %e"%lambdap)
    return B, freq, puls, rlc, champE, gammamax, rlarmor, lambdap


def rad_to_hms(angle_radians):
    '''
       Routine pour convertir des radians en h:min:sec
        angle_radians est un angle en radians > 0
        h,m,s = Rad_to_hms(angle_radians)
    '''
    angle = abs(angle_radians)
    heure = pi / 12.
    minute = heure / 60.
    seconde = minute / 60.
    h = int(angle / heure)
    m = int( (angle - h * heure) / minute)
    s = ( angle - h*heure - m * minute ) / seconde
    return h,m,s


def hms_to_rad(h, m, s):
    '''
       Routine pour convertir des h:min:sec en radians 
        h, m, s are resp. hours minutes and seconds > 0
        angle_radians = hsms_rad(h,m,s)
    '''
    heure = pi / 12.
    minute = heure / 60.
    seconde = minute / 60.
    return h*heure + m * minute + s * seconde


def rad_to_dms(angle_radians):
    '''
       Routine pour convertir des radians en degré:minute de degre:second de degré
        angle_radians est un angle en radians > 0
        d,m,s = Rad_to_hms(angle_radians)
    '''
    angle = abs(angle_radians)
    degre = pi / 180.
    minute = degre / 60.
    seconde = minute / 60.
    d = np.int(np.floor(angle / degre))
    m = np.int(np.floor( (angle - d * degre) / minute))
    s = ( angle - d*degre - m * minute ) / seconde
    return d,m,s

def dms_to_rad(d,m,s):
    '''
       Routine pour convertir des h:min:sec en radians 
        d, m, s are resp. degrees minutes of degree and seconds of degree > 0
        angle_radians = dms_rad(d,m,s)
    '''
    degre = pi / 180.
    minute = degre / 60.
    seconde = minute / 60.
    return d*degre + m * minute + s * seconde





def Position_with_proper_motion(ra, dec, ra1, dec1, delta_t) :
    '''
        ra1 in mas of degree  / yr
        dec1 in mas of degree / yr
        delta_t in MJD
    '''
    # Conversion from mas/yr to radians of arc per day */
    convert = 1.0/1000.0/60.0/60.0 * pi/180.0 /365.25;

    r0 = np.array([np.cos(ra) * np.cos(dec), np.sin(ra) * np.cos(dec), np.sin(dec)])
    alpha = np.array([-np.sin(ra), np.cos(ra),0.])
    delta = np.array([-np.cos(ra) * np.sin(dec), -np.sin(ra) * np.sin(dec), np.cos(dec)])
    rfinal = r0 + ( alpha * ra1 + delta * dec1) * delta_t * convert
    dirfinal = rfinal / np.sqrt((rfinal**2).sum())
    cdec = np.sqrt(dirfinal[0]**2 + dirfinal[1]**2)
    sdec = dirfinal[2]
    decfinal = inversetrigo(cdec,sdec)

    cra = dirfinal[0] / cdec
    sra = dirfinal[1] / cdec
    rafinal = inversetrigo(cra, sra)

    return rafinal, decfinal
    
    


def Parallax_to_distance(px, unit='ly'):
    '''
        px is an annual parallax given in in mas
        if unit is 'ly', give the corresponding distance in light years, otherwise is parsecs
    ''' 
    distpc = 1000./ px 
    if unit =='ly' :
        return distpc * meterpc / meterly
    else :
        return distpc

