# SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later

#coding:utf8
import numpy as np

pi = np.arccos(-1.)

def rotation_2D(vector, thet):
    '''y
        ^  _
        |  /| x'
        | /
        |/_______> x
        thet = (x,x')
    vector is expressed in the x' base 
    out is expressed in the x base
    example : for thet = pi/2 : x' becomes y : (1,0) -> (0,1)
    '''
    out = np.zeros(2)
    costhet = np.cos(thet)
    sinthet = np.sin(thet)
    out[0] = vector[0] * costhet - vector[1] * sinthet
    out[1] = vector[0] * sinthet + vector[1] * costhet
    return out 



def inversetrigo( cosv , sinv) :
    '''
        Return angle v with cosine cosv and sine sinv. Adapted from c++ version in utilities.cpp (project triple system)
    '''

    v = 0. 
    if ( abs(sinv ) > 1. or abs(cosv) > 1. ) :
        print("inversetrigo : Sinus ou cosinus supérieur à 1 !!!sin = %e et cos= %e "%( sinv ,cosv))  
        sinv = min(sinv, 1.)
        sinv = max(sinv, -1.)
        cosv = min(cosv, 1.)
        cosv = max(cosv, -1.)
        print("inversetrigo : utilise sin = %e et cos = %e"%(sinv, cosv)) 
        
    if (sinv >= 0. and cosv >= 0. ) :
        v = np.arcsin( sinv ) 
    elif (sinv >= 0. and cosv < 0. ) :
        v = np.arccos( cosv ) 
    elif (sinv  < 0. and cosv >= 0. ) :
        v = np.arcsin( sinv ) 
    elif (sinv < 0. and cosv < 0. ) :
        v = pi + np.arctan( sinv / cosv ) 
    return v 


def sphere2cart(coord_sphere):
    '''
        convert spherical coordinates into cartesian coordinates assuming thet is taken with respect to z with :
        coord_sphere = array(r,thet,phi)
        coord_cart = array(x,y,z)
    '''
    sinthet = np.sin(coord_sphere[1])
    return coord_sphere[0] * np.array([sinthet * np.cos(coord_sphere[2]), sinthet * np.sin(coord_sphere[2]), np.cos(coord_sphere[1])])


def sphere2cart_vector(vector_sphere, coord_sphere):
    '''
        convert spherical coordinates components of a vector into cartesian coordinates components at location coord_sphere assuming thet is taken with respect to z with :
        coord_sphere = array(r,thet,phi)
        coord_cart = array(x,y,z)
        return vector_cart 
    '''
    vector_cart= np.zeros(3)

    costhet = np.cos(coord_sphere[1])
    sinthet = np.sin(coord_sphere[1])
    cosphi = np.cos(coord_sphere[2])
    sinphi = np.sin(coord_sphere[2])

    vector_cart[0] = (vector_sphere[0] * sinthet * cosphi +
                      vector_sphere[1] * costhet * cosphi -
                     vector_sphere[2] * sinphi )
    vector_cart[1] = ( vector_sphere[0] * sinthet * sinphi +
                     vector_sphere[1] * costhet * sinphi +
                     vector_sphere[2] * cosphi )
    vector_cart[2] = ( vector_sphere[0] * costhet -
                     vector_sphere[1] * sinthet )

    return vector_cart


def cart2sphere(coord_cart) :
    '''
        convert cartesian coordinates into spherical coordinates assuming thet is taken with respect to z with :
        coord_cart = array(x,y,z)
    '''
    eps = 0.0000000000001
    r = np.sqrt( (coord_cart**2).sum() )
    if r > eps :
        costhet = coord_cart[2] / r
        thet = np.arccos(costhet)
        sinthet = np.sqrt( 1. - costhet**2)
        if sinthet > eps :        
          #  rcyl = r* sinthet 
            tanphi = coord_cart[1]/coord_cart[0] 
            if abs(tanphi) > 1./eps :
                phi = 0.5*pi *np.sign(coord_cart[1])
            elif abs(tanphi) < eps :
                phi = (1. + np.sign(coord_cart[0]) ) * pi
            else :
                cosphi = np.sign(coord_cart[0])/(np.sqrt(1. + tanphi**2))
                sinphi = np.sign(coord_cart[1])/(np.sqrt(1. + 1./tanphi**2))
                
                #cosphi = coord_cart[0] / rcyl
                #sinphi = coord_cart[1] / rcyl
                phi = inversetrigo(cosphi, sinphi)
        else :
            print(" warning in cart2sphere : sinthet too small !") 
            phi = 0.
    else :
        print("warning in cart2sphere : r too small !") 
        thet = 0.
        phi = 0.
    return np.array([r, thet, phi])



def dichotomie(func, xmin, xmax, args=(), goal = 0., maxit =100, epsrelmax = 0.00001, epsx = -1, verbose = 0):
    '''
        Résout par dichotomie func(x,*args) = goal et renvoie le résultat pour x. 
        xmin et xmax : points de départ
        maxit  : nb max d'itérations
        epsmax : précision à atteindre
        epsx   : if > 0.,  altenative stopping criterion. Stop if the distance between the approximate result xx and the exact one is less than epsx 
        verbose : > 0 then say stuff
    '''
    xxmin = xmin 
    xxmax = xmax
    funcxxmax = func(xxmax, *args)
    funcxxmin = func(xxmin, *args) 
    sens = 0
    if funcxxmin > goal and funcxxmax < goal : 
        sens = 1
    elif funcxxmin < goal and funcxxmax > goal :
        sens = -1
    else : 
        print(" Wrong range ! ")
        return xmin

    xx = (xxmax + xxmin) / 2.
    funcxx = func(xx, *args) 
    i = 0
    
    if goal == 0. :
        epsmax = epsrelmax
    else :
        epsmax = epsrelmax * abs(goal)

    while(abs(funcxx - goal) > epsmax and i < maxit) :
        if (funcxx - goal ) > 0. :
            if sens == 1 :
                xxmin = xx
                funcxxmin = funcxx 
            else :
                xxmax = xx
                funcxxmax = funcxx
        else :
            if sens == 1 :
                xxmax = xx 
                funcxxmax = funcxx
            else :
                xxmin = xx
                funcxxmin = funcxx 
        if verbose > 0 :
            print("dichotomie state iteration, xmin, xmax", i, xxmin, xxmax)
        #funcxxmax = func(xxmax,*args) 
        #funcxxmin = func(xxmin, *args) 
        xx = (xxmax + xxmin) / 2.
        funcxx = func(xx, *args)

        if abs(xxmin -xxmax) < epsx : # alternative stopping criterion 
            epsmax = 2.* abs(funcxx - goal)

        i += 1

    if i == maxit : print("Warning : dichotomie reached maxit !", xmin, xmax, goal, abs(funcxx - goal))

    return xx


def Print_value_with_error_bar(value,error) :
    '''
        Prints in Latex format value as x.xxx(yy)\cdot 10^(z) where yy are two digits 
        corresponding to error and x and z are the base and exponent of value.
    '''
    expo = np.int(np.floor(np.log10(np.abs(value))))
    nf = np.int(np.floor(np.log10(error)))
    base = np.sign(value) * value/ 10**expo
    str = "${:.{prec}f}".format(base, prec= -(nf - expo) -1)
    str+="({:d})".format(np.int(error *10 /10**nf))
    str+= r"\cdot 10^{" + "{:d}".format(expo) + r"}$"
    print(str)
    return


def Print_numbers(number, chiffre_significatif=2, limite_notasci =0, epsint = 10.**(-14) ):
    '''
        Print number with chiffre_significatif significant numbers and scientific notation
        limite_notasci gives the number of numbers needed to print under which a floating point notation is prefered
    '''
    if number == 0 :
        return '{num:{width:d}f}'.format(num = number, width = chiffre_significatif)
    lognb = np.log10(number)
    if abs(lognb ) >= limite_notasci :
        if lognb > 0 :
            power = np.int(np.floor(lognb) )
        else :
            power = -np.int(np.floor(abs(lognb)))  - 1
        sub10part = number / 10.**(power) 
        if abs(sub10part - 10) < epsint :
                power += 1
                sub10part =1.
        if power == 0 :
            return '{num:.{prec:d}f}'.format(num=sub10part,  prec = max(chiffre_significatif-1,0))
        else :
            if chiffre_significatif == 0 :
                return "10^{{{expo:d}}}".format(expo = power)
            else :
                return "{num:.{prec:d}f}".format(num=sub10part,  prec = chiffre_significatif-1,) +   '\\' + "cdot 10^{{{expo:d}}}".format(expo = power)
    else :
        return '{num:{width:d}f}'.format(num = number, width = chiffre_significatif)
