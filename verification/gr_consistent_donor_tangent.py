"""Counterexample candidate: differentiate species on the selected donor too."""
from pathlib import Path
from types import FunctionType

import numpy as np
import gr_upwind_iteration_repair as core

m,e,ld=core.m,core.e,core.ld
OUT=m.BASE/'step42-consistent-donor'


def selected_species(star,candidate,previous,older,h,coefficients,metric):
    _,c1,c2=coefficients
    v=candidate[:,2]
    B=metric['a']*star.reference['rho']*np.exp(candidate[:,0])/np.sqrt(1-v*v)
    velocity=star.faces(metric['N']*v/metric['a'],odd=True)
    donor=star.selected_donor
    flux=velocity*m.prior.upwind(B,donor)
    effective=-c1*previous[1]['B']-c2*older[1]['B']
    assert np.all(effective>0)
    left=h*e.C*4*np.pi*star.rf[:-1]**2*np.where(donor[:-1]>=0,flux[:-1],0)/star.volume
    right=h*e.C*4*np.pi*star.rf[1:]**2*np.where(donor[1:]<0,-flux[1:],0)/star.volume
    band=np.zeros((3,star.n))
    band[1]=1+(left+right)/effective
    band[0,1:],band[2,:-1]=-right[:-1]/effective[:-1],-left[1:]/effective[1:]
    changes=np.diff(star.base[:,5:],axis=0)
    source=left[:,None]*np.vstack([np.zeros(26),-changes])+right[:,None]*np.vstack([changes,np.zeros(26)])
    source+=-c1*previous[1]['B'][:,None]*previous[0][:,5:]-c2*older[1]['B'][:,None]*older[0][:,5:]
    return m.solve_banded((1,1),band,np.asarray(source/effective[:,None],float)).astype(ld)


evaluate=FunctionType(m.method.CompositionTangent.evaluate.__code__,dict(vars(m.method),species=selected_species))


class ConsistentDonorTangent(core.PredictedDonorTangent):
    def evaluate(self,delta):
        previous,older,h,coefficients=self.composition_context
        self.projected_anchor=selected_species(self,self.anchor,previous,older,h,coefficients,self.linearization)
        return evaluate(self,delta)


context=dict(vars(core),OUT=OUT,__file__=__file__,PredictedDonorTangent=ConsistentDonorTangent)
correction=FunctionType(core.correction.__code__,context)
context['correction']=correction
finish=FunctionType(core.finish.__code__,context)
context['finish']=finish
check=FunctionType(core.check.__code__,context)


if __name__=='__main__':
    check()
