import numpy as np
import json
from pathlib import Path
import def_gr_matrix_propagation as b
o=b.OUT
sigma=12.;beta=512.;count=4096
theta=2*np.pi*(np.arange(count)+.5)/count
z=np.exp(1j*theta);s=sigma+beta*(1+z)/(1-z)
times=b.prior.TIMES
L=np.zeros((len(times),2048));L[:,0]=np.exp(-beta*times);L[:,1]=(1-2*beta*times)*L[:,0]
for n in range(1,2047):L[:,n+1]=((2*n+1-2*beta*times)*L[:,n]-n*L[:,n-1])/(n+1)
L*=np.exp(sigma*times[:,None])
def inverse(fp,n):
 G=2*beta/(1-z)*fp
 a=np.fft.fft(G)/count*np.exp(-1j*np.pi*np.arange(count)/count)
 return L[:,:n]@a[:n].real
errors={}
for omega in [3.,50.,300.,600.,800.,2000.]:
 fp=1/(s*(s*s+omega*omega));exact=(1-np.cos(omega*times))/omega**2
 errors[str(omega)]=[float(np.max(abs(inverse(fp,n)[1:]-exact[1:]))*omega**2) for n in [512,1024,2048]]
result=dict(classification='Counterexample candidate',method='Weeks Laguerre expansion, no nonlinear Q-D table or reduced-matrix spectrum',sigma=sigma,beta=beta,contour_samples=count,coefficients=[512,1024,2048],errors=errors,new_GR_solves=0)
b.write(o/'weeks-control.json',result);print(json.dumps(result),flush=True)
