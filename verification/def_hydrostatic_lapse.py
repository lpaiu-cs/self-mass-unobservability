"""Recover the lapse on the accepted mechanical solution, without refitting."""
from pathlib import Path
import json
import time
import numpy as np
import def_hydrostatic_background as h

OUT=h.OUT/'absolute-shoot'


class Structure(h.Structure):
    def radial(self,y,i):
        base=super().radial(y[:5],i);r=y[0]*self.R;m=y[1]*self.B
        phi=self.phi0*(1+self.mu*y[3]);v=self.phi0*self.mu*y[4]/self.R
        p,_,_,_=self.state(y[2],i);pE=np.exp(2*self.beta*phi*phi)*p
        nr=m/(r*r*(1-2*m/r))+4*np.pi*r*pE/(1-2*m/r)+r*v*v/2
        return np.r_[base,nr*base[0]*self.R]


def first_integral(s,y,i):
    _,_,_,v=s.state(y[2],i);p=np.exp(y[2]);rho=np.exp(v[0]);cx=s.data['CX'][i]
    phi=s.phi0*(1+s.mu*y[3])
    return y[5]+s.beta*phi*phi/2+np.log1p((v[2]+p/rho/(h.gr.C*100)**2)/cx)


def run():
    assert not (OUT/'lapse.json').exists()
    audit=json.loads((OUT/'native-audit.json').read_text());assert audit['passed']
    files=[Path(__file__),OUT/'background-0.001.npz',OUT/'native-audit.json',OUT/'result.json']
    h.write(OUT/'lapse-plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in files},
        claim='Integrate log-lapse on the accepted five-variable background and check the cellwise isentropic first integral log(h_J A N).',
        budget=dict(hard_timeout_seconds=30,new_native_calls=0,new_background_fits=0),
        gates=dict(replayed_state=1e-8,first_integral_absolute=1e-12),
        scope='Finite-pressure optical boundary; formal Just normalization, no physical vacuum junction or thermal stationarity.'))
    start=time.monotonic();saved=np.load(OUT/'background-0.001.npz');x=saved['parameters']
    base=h.Structure(.001);error,inner,outer,ys=base.branches(x,record=True);assert abs(error).max()<1e-8
    s=Structure(.001);n=len(s.f);nu=np.zeros(n+1);midnu=np.zeros(n)
    R=ys[0]*s.R;mu=ys[1]*s.B/R;q=.001*s.mu*x[4]
    h.mp.mp.dps=60;ext=h.exterior.exact(h.mp.mpf(mu),h.mp.mpf(q))
    nu[0]=float(h.mp.log1p(-2*h.mp.mpf(mu))/2-h.mp.mpf(q)**2*ext[3])
    defects=[];replay=[]
    for i in range(n):
        if i<s.split:
            coordinate=outer[i,0];lo=nu[i]
            if i==0:lo-=s.radial(np.r_[ys,nu[0]],0)[5]*s.B*coordinate
            initial=np.r_[outer[i,1:],lo]
            end=s.step(np.log(coordinate),np.log(s.outer[i+1]),initial,i,s.B,True)
            mid=s.step(np.log(coordinate),np.log((s.outer[i]+s.outer[i+1])/2),initial,i,s.B,True)
            nu[i+1]=end[5];midnu[i]=mid[5]
            replay.append(float(abs(end[:5]-saved['faces'][i+1]).max()))
            defects.append(abs(first_integral(s,end,i)-first_integral(s,initial,i)))
        else:
            j=n-1-i;coordinate=inner[j,0];initial=np.r_[inner[j,1:],0.]
            end=s.step(np.log(coordinate),np.log(s.inner[i]),initial,i,s.B,False)
            mid=s.step(np.log(coordinate),np.log((s.inner[i]+s.inner[i+1])/2),initial,i,s.B,False)
            lo=nu[i]-end[5];nu[i+1]=lo;midnu[i]=lo+mid[5]
            if i==n-1:
                # Regular centre correction: nu(r0)-nu(0), through r0^2.
                p,en,_,_=s.state(x[0],i);phi=.001*(1+s.mu*x[3]);A4=np.exp(2*s.beta*phi**2)
                nu[-1]-=2*np.pi/3*A4*(en+3*p)
            replay.append(float(abs(end[:5]-saved['faces'][i]).max()))
            defects.append(abs(first_integral(s,end,i)-first_integral(s,initial,i)))
    passed=max(replay)<1e-8 and max(defects)<1e-12
    np.savez_compressed(OUT/'lapse.npz',nu_faces=nu,nu_mid=midnu,cell_first_integral_defects=defects)
    value=dict(classification='Counterexample candidate',passed=bool(passed),
        max_replayed_state_difference=max(replay),max_cell_first_integral_defect=max(defects),
        centre_log_lapse=float(nu[-1]),surface_log_lapse=float(nu[0]),seconds=time.monotonic()-start,
        native_calls=0,new_background_fits=0,physical_atmosphere=False,thermal_stationarity=False)
    h.write(OUT/'lapse.json',value);print(json.dumps(value),flush=True);assert passed


if __name__=='__main__':run()
