# SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later

# Guillaume VOISIN , LUTh, Observatoire de Paris - PSL (guillaume.voisin@obspm.fr astro.guillaume.voisin@google.com)  

import numpy as np
from matplotlib import pyplot as pplot

class Read_and_Sieve:
    def __init__(self, filename):
        '''
        Sieve and trim mcmc resut file according to lnposterior conditions. 
        
        filename : path to the mcmc result file
        '''
        self.lnpost = []
        self.lines = []
        file = open(filename,"r")
        line = file.readline()
        while (line != '') :
            sline = line.split()
            if (sline[0][0] !='#') :
                self.lines.append(line)
                self.lnpost.append(np.float64(sline[-1]))
            line = file.readline()
        file.close()

        ##

        return

    def Sieve(self, threshold):
        '''
        Remove from chain all elements with lnposterior < threshold 
        The sieved chain is stored into self.lnpost_sieved for the lnposteriors 
        and self.lines_sieved for the unparsed file lines. 
        '''
        self.lnpost_sieved = []
        self.lines_sieved = []
        for i in range(len(self.lnpost)):
            if self.lnpost[i] > threshold :
                self.lnpost_sieved.append(self.lnpost[i])
                self.lines_sieved.append(self.lines[i])
        return

    def Trim_copies(self, eps = 1.e-5):
        '''
        Need to run Sieve first !
        
        Filters the chain such that any two elements have abs(lnposterior1 - lnposterior2) > eps
        
        lnposteriors are saved to self.lnpost_sieved_trimmed
        unparsed file lines are saved to self.lines_sieved
        '''
        self.lnpost_sieved_trimmed = []
        self.lines_sieved_trimmed = []
        rejected = 0
        for i in range(len(self.lnpost_sieved)):
            reject = 0
            j=0
            while (j < len(self.lnpost_sieved_trimmed) and reject == 0):
                if abs(self.lnpost_sieved_trimmed[j] - self.lnpost_sieved[i]) < eps :
                    reject = 1
                j += 1
            if reject == 0:
                self.lnpost_sieved_trimmed.append(self.lnpost_sieved[i])
                self.lines_sieved_trimmed.append(self.lines_sieved[i])
            else :
                rejected +=1
        print("Trimmed {:} copies.".format(rejected))
        return

    def Do_all(self,filename, threshold, eps=1.e-5, nlast=-1):
        '''
        Sieve then trim then save to file "filename"
        
        threshold : minimum ln_posterior = -chi2/2 to conserve (see self.Sieve)
        eps : minimum distance between any two elements kept after trimming (see self.Trim_copies)
        nlast :  if nlast > 0, save only the last "nlast" lines. Otherwise save all. (see self.Save_sieved_trimmed)
        '''
        self.Sieve(threshold)
        self.Trim_copies(eps=eps)
        self.Save_sieved_trimmed(filename, nlast=nlast)
        return

    def Plot_chi2(self, dof = -1, bins=50, start = 0, end = -1):
        '''
        Plot the histogram of chi2 = -2*self.lnposterior

        dof : number of degrees of freedom to overlay the ideal chi2 distribution. If < 0, then nothing is plotted.
        bins : histogram bins
        start : rank of the first chi2 to plot
        end : rank of the last chi2 2 plot. If < 0, is treated like a numpy array.
        '''
        fig = pplot.figure()
        plt = fig.add_subplot(111)

        if end < 0 :
            eend = len(self.lnpost) - end
        else :
            eend = end

        chi2counts, chi2bins, patches = plt.hist(abs(np.array(self.lnpost[start:eend])*2), bins=bins)

        if dof > 0:
            #dof = self.res.shape[1] #ndatapoints - self.res.shape[1] # number of degrees of freedom
            chi2bins = 0.5*(chi2bins[:chi2bins.size-1] + chi2bins[1:])
            chi2bins_mean = np.average(chi2bins, weights=chi2counts)
            chi2bins_mean = min(chi2bins_mean, chi2bins[0]+dof-0.1)
            chi2s = chi2bins - chi2bins_mean + dof #*temperature
            chi2distro = np.log(chi2s)*(dof/2. -1.) + (-chi2s/2.)
            chi2distro = np.exp(chi2distro)
            chi2distro *= chi2counts.max() / chi2distro.max() # rough scaling
            plt.plot(chi2bins, chi2distro, label=r'Estimate of $\chi^2$ distribution {:d} dof'.format(dof))
            fig.legend(loc=1)

        plt.set_xlabel(r'$\chi^2$')
        plt.set_ylabel('Counts')
        plt.ticklabel_format(style='sci', scilimits=(0,0), axis='y')

        fig.tight_layout(rect=(0,0,1,0.90), pad =1.11)
        fig.set_tight_layout(True)

        fig.show()

        self.fig = fig
        self.plt = plt
        return

    def Save_sieved(self, filename, nlast = -1):
        '''
        Save the sieved chain to filename
        if nlast > 0, save only the last "nlast" lines. Otherwise save all. 
        
        Warning : the file header needs to be added by hand. 
        '''
        file = open(filename,"w")
        nlines = len(self.lines_sieved)
        if nlast >0 :
            start = nlines - nlast
        else :
            start =0
        for i in range(start, nlines): #for line in self.lines_sieved:
            line = self.lines_sieved[i]
            file.write(line)
        file.close()
        print("Warning : the file header needs to be added by hand. ")
        return

    def Save_sieved_trimmed(self, filename, nlast=-1):
        '''
        Save the sieved chain to filename
        if nlast > 0, save only the last "nlast" lines. Otherwise save all. 
        
        Warning : the file header needs to be added by hand. 
        '''
        file = open(filename,"w")
        nlines = len(self.lines_sieved_trimmed)
        if nlast >0 :
            start = nlines - nlast
        else :
            start =0
        for i in range(start, nlines):# in self.lines_sieved_trimmed:
            line = self.lines_sieved_trimmed[i]
            file.write(line)
        file.close()
        print("Warning : the file header needs to be added by hand. ")
        return
