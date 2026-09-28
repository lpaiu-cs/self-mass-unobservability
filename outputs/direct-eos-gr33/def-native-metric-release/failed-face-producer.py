"""Prepare a conserved GR source with a common physical inner-face state.

Counterexample candidate: repair the reconstructed background at the gas/bulk
interface before interpreting its numerical boundary transfer as a GR source.
The physical gravity and photon source terms are retained.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect
import json
import signal
import time
import numpy as np
import sympy as s
import def_native_energy_flow as old

OUT=old.OUT.parent/'def-native-metric-release'
C=old.C
write=old.write


def symbolic():
    r,m,A,alpha,E,P,phi,f,fr,J,dE,dP=s.symbols('r m A alpha E P Phi f fr J dE dP',nonzero=True)
    b=1-2*m/r;mr=4*s.pi*r*r*A**4*E+r*r*b*phi**2/2
    nu=m/(r*r*b)+4*s.pi*r*A**4*P/b+r*phi**2/2
    lam=(mr/r-m/r**2)/b
    phir=-(2/r+nu-lam)*phi+4*s.pi*alpha*A**4*(E-3*P)/b
    # All background/perturbation stresses here are in geometric units.
    dm=r*r*b*phi*f+J
    Jr=4*s.pi*r*r*A**4*(dE+(E+P)*(3*alpha+r*phi)*f)-r*phi**2*J
    dmr=sum(s.diff(dm,x)*dx for x,dx in [(r,1),(m,mr),(phi,phir),(f,fr),(J,Jr)])
    target=4*s.pi*r*r*A**4*(dE+4*alpha*E*f)-r*phi**2*dm+r*r*b*phi*fr
    assert s.simplify(dmr-target)==0
    dl=dm/(r*b);dlr=sum(s.diff(dl,x)*dx for x,dx in [(r,1),(m,mr),(phi,phir),(f,fr),(J,Jr)])
    dnr=dm/(r*r*b*b)+r*phi*fr+4*s.pi*r*A**4/b*(dP+4*alpha*P*f+2*P*dl)
    expected=(2/(r*r*b*b)+8*s.pi*A**4*(P-E)/b**2)*dm+4*s.pi*r*A**4/b*(dP-dE)+16*s.pi*alpha*r*A**4/b*(P-E)*f
    assert s.simplify(dnr-dlr-expected)==0
    assert s.simplify(r*phi**2-(nu+lam)+4*s.pi*r*A**4*(E+P)/b)==0
    return dict(classification='Proven',passed=True,
        mass='dm=r^2*b*Phi*f+J; J_prime+r*Phi^2*J=4*pi*r^2*A^4*[dE+(E+P)*(3*alpha+r*Phi)*f].',
        metric_gradient='dnu_prime-dlambda_prime=[2/(r^2*b^2)+8*pi*A^4*(P-E)/b^2]*dm+4*pi*r*A^4*(dPr-dE)/b+16*pi*alpha*r*A^4*(P-E)*f/b.',
        integrating_factor='If mu_prime/mu=r*Phi^2, then (ln(mu*sqrt(b)/N))_prime=-4*pi*r*A^4*(E+P)/b. The factor is constant only in vacuum.',
        scope='Linearized polar-areal scalar and Einstein constraints about the declared static isotropic background. These identities do not evolve matter, prove the radiation background stationary, or justify replacing a numerical interface flux by a physical one.')


namespace=dict(vars(old),OUT=OUT)
rhs_source=inspect.getsource(old.Flow.rhs)
before="""        ghost=np.array(m.left);ghost[1]=np.interp(t,m.hist_t,m.hist_v)
        ext=np.column_stack([ghost,V,np.zeros(3)]);slope=np.zeros_like(ext)
        slope[:,1:-1]=bank.task.minmod(ext[:,1:-1]-ext[:,:-2],ext[:,2:]-ext[:,1:-1])
        L=ext[:,:-1]+slope[:,:-1]/2;R=ext[:,1:]-slope[:,1:]/2"""
assert rhs_source.count(before)==1
rhs_source=rhs_source.replace(before,'        L,R=self.reconstruct(V,t)')
import textwrap
exec(compile(textwrap.dedent(rhs_source),__file__,'exec'),namespace)


class Flow(old.Flow):
    rhs=namespace['rhs']
    run=FunctionType(old.Flow.run.__code__,namespace)
    primitive=FunctionType(old.Flow.primitive.__code__,namespace)

    def __init__(self,n):
        super().__init__(n);m=self.base
        self.eos=old.temperature.parent.EOS(old.OUT/'runtime-columns.npz');self.eos.fan=SimpleNamespace(calls=0)
        self.background_cell=np.array([self.initial[0],np.zeros(n),self.seed.copy()])
        env=m.env;d=m.eos.d;r=m.R+m.xf/m.As
        density=np.exp(np.interp(np.minimum(r,m.R),env['r'],np.log(env['rho'])))/m.eos.rho0
        entropy=np.interp(np.minimum(r,m.R),d['initial_r'],d['initial_sigma'])
        theta=np.log(m.eos(density,entropy)[3]);vacuum_theta=float(np.log(m.eos(np.array([0.]),np.array([0.]))[3][0]))
        self.background_left=np.array([density,np.zeros(n+1),theta]);self.background_right=self.background_left.copy()
        self.background_left[:,m.xf>0]=np.array([0.,0.,vacuum_theta])[:,None]
        self.background_right[:,m.xf>=0]=np.array([0.,0.,vacuum_theta])[:,None]

    def reconstruct(self,V,t):
        m=self.base;delta=V-self.background_cell;ghost=np.array([0.,np.interp(t,m.hist_t,m.hist_v),0.])
        ext=np.column_stack([ghost,delta,np.zeros(3)]);slope=np.zeros_like(ext)
        slope[:,1:-1]=old.bank.task.minmod(ext[:,1:-1]-ext[:,:-2],ext[:,2:]-ext[:,1:-1])
        L=self.background_left+ext[:,:-1]+slope[:,:-1]/2
        R=self.background_right+ext[:,1:]-slope[:,1:]/2
        assert min(L[0])>=-1e-13 and min(R[0])>=-1e-13,'Positive reconstructed gas'
        L[0]=np.maximum(L[0],0);R[0]=np.maximum(R[0],0)
        return L,R


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='a0db40a17',
        claim='Before applying conserved matter to Einstein constraints, test and repair the inner-face background reconstruction without removing physical gravity or radiation forces. The pressure release at the old surface remains a real gas/vacuum jump.',
        repair='Reconstruct perturbations around the saved nonuniform gas background. Both sides of each interior face use its common native background rho,T; the old surface retains gas on its left and vacuum on its right. Keep all geometric and scattering source terms.',
        paths=[896,1792],horizon_seconds=.0034344311179287023,domain_m=[-200,1200],native_calls=0,
        seconds=100,CPU_threads=1,memory_GB=2,
        gates=dict(energy_ledger=1e-8,baryon_ledger=1e-10,outside_mass_refinement=.02,trace_refinement=.02),
        stop='Use only the frozen final native EOS support; stop on a domain, conservation, refinement or budget failure. No new native bank, resolution or time extension.',
        bindings={str(p.relative_to(old.bank.old.ROOT)):old.bank.old.photons.digest(p) for p in [Path(__file__),Path(old.__file__),old.OUT/'runtime-columns.npz',old.OUT/'cells-1792.npz']}))
    write(OUT/'symbolic.json',symbolic());(OUT/'reused-rhs.py').write_text(rhs_source)
    original=old.Flow(1792);original.eos=old.temperature.parent.EOS(old.OUT/'runtime-columns.npz');original.eos.fan=SimpleNamespace(calls=0)
    fixed=Flow(1792);rows=[]
    for label,m in [('original',original),('common_face_background',fixed)]:
        rate,ledger,dt=m.rhs(m.initial,0.);rho=m.initial[0];p,u,*_=m.eos(rho,m.seed);h=m.eos.cx+u+p/np.maximum(rho,m.eos.floor)
        acceleration=rate[1]/np.maximum(rho*h,1e-100)*C;inside=m.base.x<0
        rows.append(dict(method=label,first_cells_acceleration_cm_s2=acceleration[:8].tolist(),
            interior_acceleration_median_cm_s2=float(np.median(acceleration[4:100])),initial_baryon_flux=float(ledger[0]),initial_energy_flux=float(ledger[1])))
    write(OUT/'initial-face-control.json',dict(classification='Counterexample candidate',rows=rows,physical_sources_subtracted=False))
    print(json.dumps(rows),flush=True)


def run():
    assert not (OUT/'result.json').exists();start=time.monotonic();signal.alarm(100)
    spec=json.loads((OUT/'plan.json').read_text())
    for p,h in spec['bindings'].items():assert old.bank.old.photons.digest(old.bank.old.ROOT/p)==h,p
    coarse=Flow(896).run('cells-896');elapsed=time.monotonic()-start
    forecast=5.2*coarse['seconds']+5.5
    write(OUT/'measured-budget.json',dict(classification='Counterexample candidate',coarse_seconds=coarse['seconds'],elapsed_seconds=elapsed,
        forecast_last_path_seconds=forecast,remaining_seconds=100-elapsed,assumption='4x cell-squared cost with30percent margin and5.5seconds setup, using the just-measured common-face producer. Fine-grid speed remains unmeasured.'))
    assert elapsed+forecast<100,'Remaining path exceeds budget'
    fine=Flow(1792).run('cells-1792')
    errors={k:abs(coarse[k]/fine[k]-1) for k in ['gas_outside_original_radius_g','integrated_trace_energy_erg']}
    write(OUT/'result.json',dict(classification='Counterexample candidate',passed=max(errors.values())<.02,
        refinement=errors,seconds=time.monotonic()-start,common_background_face_applied=True,physical_source_subtracted=False,
        metric_constraints_derived=True,full_GR_scalar_feedback=False,final_charge_solved=False,full_goal_complete=False))
    signal.alarm(0);print((OUT/'result.json').read_text(),flush=True)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
