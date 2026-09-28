"""Counterexample candidate: close the temperature-only conductive GR loop.

The native EOS tangent and GR response determine Eulerian temperature. Only
that temperature is returned to the original pole drive. Transport geometry,
conductivities and rates stay frozen; this is not the full nonlinear channel.
"""
from pathlib import Path
import argparse
import json
import resource
import signal
import time
import numpy as np
from scipy.sparse import coo_matrix, diags, csc_matrix
import def_gr_interface_patch as patch
import def_gr_cached_resolvent as cached

go=patch.go
OUT=patch.OUT.parent/'def-gr-temperature-feedback'
write=go.write


class Problem(cached.Problem):
    def __init__(self,degree=4):
        super().__init__(degree)
        m=self.model; heat=m.heat; d=heat.d; r=heat.rnative
        assert np.array_equal(r,m.original.native)
        V,D=m.evaluation(r); point=m.bg.sample(r)
        N,a=heat.geometry.metric(r); mass,p,phi,v=[point[k] for k in ['m','p','phi','v']]
        b=1-2*mass/r; A4=np.exp(-8*phi**2); alpha=-4*phi
        ad=point['adiabatic_T_rho']
        dlz=r*r*v*v/2-4*np.pi*r*r*A4*p/b-mass/(r*b)
        # Same piecewise-linear background lnT as the native thermal input.
        ids=np.clip(np.searchsorted(r,r,side='right')-1,0,len(r)-2)
        slope=np.diff(d['lnT'][::-1])[ids]/np.diff(r)[ids]
        Q=-diags(ad*r)@D[0]-diags(ad*(3+dlz+3*alpha*r*v)+r*slope)@V[0]
        Q-=diags(ad*(r*v+3*alpha))@V[1]
        _,loss,J=go.task.fem.source_points(heat,point)
        gr=go.task.fem.base.task.h.gr; geo=gr.G*.1*heat.geometry.R**2/gr.C**4
        rho=d['raw'][::-1,0]; cvT=d['thermo'][::-1,5]
        TE=-diags(ad/(r*b))@J-diags(1/(rho*geo*cvT))@loss
        theta=(d['A']*d['N']*np.exp(d['lnT']))[::-1]
        self.Tq=Q.astype(np.clongdouble); self.TE=TE.astype(np.clongdouble)
        self.theta=theta; self.V=V; self.D=D; self.point=point
        bank=np.load(go.task.BANK/'fine-bank.npz'); faces=bank['faces']; n=len(r)
        diff=coo_matrix((np.tile([-1.,1.],len(faces)),
             (np.repeat(np.arange(len(faces)),2),np.c_[n-faces,n-1-faces].ravel())),shape=(len(faces),n)).tocsr()
        Nf=np.sqrt(d['N'][faces-1]*d['N'][faces]); Af=np.sqrt(d['A'][faces-1]*d['A'][faces])
        af=np.sqrt(d['metric'][faces-1]*d['metric'][faces])
        pref=-4*np.pi*d['faces_cm'][faces]**2*Nf*Af**2/af*1e5/np.diff(d['radius_cm'])[faces-1]
        self.conductance=(pref[:,None]*bank['mode_K_SI']).astype(np.longdouble)
        self.Gq=(diff@diags(theta)@Q).astype(np.clongdouble)
        self.GE=(diff@diags(theta)@(TE-Q@m.H)).astype(np.clongdouble)
        self.background_difference=np.asarray(diff@theta)
        assert np.allclose(self.conductance*self.background_difference[:,None],heat.amplitude,rtol=3e-15,atol=0)
        self.loop_residual=0.; self.iterations=0; self.second_to_first=0.; self.feedback_relative=0.

    def transform_pair(self,z):
        m=self.model; tc=np.longdouble(m.original.radiation.geometry.tc)
        gain=np.sum(self.conductance*tc*self.lam/(z*(z+self.lam)),axis=1)
        E0=np.zeros(len(m.heat.edges),np.clongdouble)
        E0[m.heat.face_ids]=np.sum(self.numerator/(z*z*(z+self.lam)),axis=1)
        data=self.Kdata+complex(z*z)*self.Mdata
        matrix=csc_matrix((data,self.indices,self.indptr),shape=self.K.shape)
        scale=np.sqrt(abs(matrix.diagonal())); inv=1/scale
        scaled=csc_matrix(((data*inv[self.indices])*inv[self.columns],self.indices,self.indptr),shape=self.K.shape)
        lu=go.splu(scaled)
        def solve(E):
            rhs=self.load@E
            u=(lu.solve(np.asarray(rhs/scale,complex))/scale).astype(np.clongdouble)
            for _ in range(3):
                defect=rhs-self.Kx@u-z*z*(self.Mx@u)
                u+=(lu.solve(np.asarray(defect/scale,complex))/scale).astype(np.clongdouble)
            defect=rhs-self.Kx@u-z*z*(self.Mx@u)
            error=float(np.max(abs(defect)/(self.absK@abs(u)+abs(z*z)*(self.absM@abs(u))+abs(rhs)+1e-100)))
            self.error=max(self.error,error); assert error<1e-9
            return u
        def feedback(u,E):
            out=np.zeros_like(E);out[m.heat.face_ids]=gain*(self.Gq@u+self.GE@E)
            return out
        u0=solve(E0); first=feedback(u0,E0); denom=max(float(abs(first).max()),1e-100)
        dE=np.zeros_like(E0); du=np.zeros_like(u0); term=first
        for iteration in range(6):
            contribution=solve(term); dE+=term; du+=contribution
            next_term=feedback(contribution,term)
            relative=float(abs(next_term).max()/denom)
            if iteration==0:self.second_to_first=max(self.second_to_first,relative)
            if relative<1e-11:break
            term=next_term
        else:raise AssertionError('Temperature feedback did not converge in six corrections')
        defect=dE-first-feedback(du,dE)
        residual=float(abs(defect).max()/denom)
        self.loop_residual=max(self.loop_residual,residual);assert residual<1e-10
        self.iterations=max(self.iterations,iteration+1)
        self.feedback_relative=max(self.feedback_relative,float(abs(dE).max()/max(abs(E0).max(),1e-100)))
        return u0-m.H@E0, du-m.H@dE, E0, dE


def control():
    import sympy as s
    z,k,c,a,b=s.symbols('z k c a b',positive=True)
    # A two-cell heat pair: Tdot=-k L T/c has no growing energy mode.
    L=s.Matrix([[1,-1],[-1,1]])
    assert s.simplify((s.Matrix([[1,1]])*L))==s.zeros(1,2)
    assert s.expand((s.Matrix([[a,b]])*L*s.Matrix([a,b]))[0]-(a-b)**2)==0
    forcing=s.Matrix([1,-1]); exact=(s.eye(2)+k*L/(c*z)).inv()*forcing
    assert s.simplify(exact-forcing/(1+2*k/(c*z)))==s.zeros(2,1)
    # u=G^-1 load E, E=E0+R(Tq*u+(TE-Tq*H)E).
    u,E,H,Tq,TE=s.symbols('u E H Tq TE')
    assert s.expand(Tq*(u-H*E)+TE*E-(Tq*u+(TE-Tq*H)*E))==0
    return dict(classification='Proven',passed=True,
        identity='Internal heat debit telescopes. Positive two-cell Fourier conduction damps the temperature difference. Temperature reconstruction uses physical q=u-H E, so the momentum lift must also be subtracted in the feedback map.',
        scope='Algebraic model identities; no uniform GR loop bound or full transport closure.')


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(90);started=time.monotonic()
    paths=[Path(__file__),Path(patch.__file__),Path(cached.__file__),Path(go.__file__),
           patch.OUT/'result.json',patch.OUT/'p4-2048.npz',go.task.BANK/'fine-bank.npz']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='ad93f3b77',
        claim='Actually return the reconstructed Eulerian temperature perturbation to all4012 existing conductive face poles and solve the resulting GR/heat response together.',
        model='Same native EOS thermal tangent and conservative heat debit; same GR mass constraint and momentum lift. delta(theta)=A0*N0*T0*delta_lnT_Eulerian, delta_lnT_Eulerian=Delta_lnT-xi*dlnT0/dr. Conductivity, proper rates, face geometry, lapse and conformal factors in the transport law remain frozen. No photon table is extended radially.',
        decision='If the solved temperature-only feedback is below1e-4 of every original native response, do not extend the time interval for this lever; prioritize missing exterior/photon channels. A larger or unresolved feedback requires a new physical error plan, not an automatic larger run.',
        numerical_scope='One p4 path only; reuse accepted baseline readouts. Measure correction directly instead of subtracting nearly equal total paths. Same4096 contour,512/1024/2048 coefficients,65 readings and0.23080495568542375s. Finite contour samples do not prove uniform contraction or continuum error.',
        gates=dict(original_linear_residual=1e-9,feedback_residual=1e-10,correction_time_absolute_to_baseline=1e-6,correction_contour_absolute_to_baseline=1e-6,baseline_endpoint_replay=1e-10,heat_balance=2e-13),
        budget=dict(pilot_seconds=90,production_seconds=600,total_seconds=720,CPU_threads=1,memory_GB=4,new_EOS_calls=0,maximum_corrections=6,automatic_expansion=False),
        bindings={str(p):go.task.digest(p) for p in paths}))
    write(OUT/'symbolic.json',control());patch.install();tick=time.monotonic();p=Problem();setup=time.monotonic()-tick
    ids=np.linspace(1,go.COUNT//2,16,dtype=int);tick=time.monotonic()
    results=[p.transform_pair(go.contour(int(k),12)[0]) for k in ids];elapsed=time.monotonic()-tick
    np.savez_compressed(OUT/'pilot.npz',ids=ids,base_q=[v[0] for v in results],correction_q=[v[1] for v in results],base_E=[v[2] for v in results],correction_E=[v[3] for v in results])
    # Use measured native long-double inversion cost, not double-BLAS speed.
    prior=json.loads((patch.PRIOR/'pilot-budget.json').read_text())
    forecast=1.4*(setup+elapsed/16*2048+prior['inversion_forecast_seconds']+30)
    endpoint=np.load(patch.OUT/'p4-2048.npz');q=endpoint['q'];E=endpoint['heat_energy']
    temp=p.Tq@q+p.TE@E
    np.savez_compressed(OUT/'temperature-endpoint.npz',radius=p.model.original.native,Eulerian_delta_lnT=temp)
    row=dict(classification='Counterexample candidate',setup_seconds=setup,transfer16_seconds=elapsed,
        forecast_seconds=forecast,seconds=time.monotonic()-started,linear_residual=p.error,feedback_residual=p.loop_residual,
        corrections=p.iterations,second_to_first=p.second_to_first,sampled_heat_feedback_relative=p.feedback_relative,
        maximum_saved_Eulerian_delta_lnT=float(abs(temp).max()),
        assumption='Sixteen representative full coupled resolvents scaled to2048 plus saved long-double native inversion,30s output allowance and40 percent margin. No uniform spectral bound.')
    write(OUT/'pilot.json',row);signal.alarm(0);print('FEEDBACK PILOT',json.dumps(row),flush=True)


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for path,h in plan['bindings'].items():assert go.task.digest(Path(path))==h,path
    pilotrow=json.loads((OUT/'pilot.json').read_text());assert pilotrow['forecast_seconds']<600
    signal.alarm(600);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)))
    started=time.monotonic();patch.install();p=Problem();m=p.model;nr=len(m.original.native)
    half=np.zeros((go.COUNT//2+1,2*nr),np.clongdouble);L=go.laguerre(12)
    ell=np.zeros(go.COUNT,np.longdouble);ell[:2048]=L[-1];kernel=go.fft(ell)[:go.COUNT//2+1]/go.COUNT
    base_end=np.zeros(m.size,np.longdouble);corr_end=base_end.copy();dEend=np.zeros(len(m.heat.edges),np.longdouble)
    pilot=np.load(OUT/'pilot.npz');saved={int(k):i for i,k in enumerate(pilot['ids'])}
    for k in range(1,go.COUNT//2+1):
        z,factor=go.contour(k,12)
        if k in saved:
            i=saved[k];q,dq,E,dE=[pilot[name][i] for name in ['base_q','correction_q','base_E','correction_E']]
        else:q,dq,E,dE=p.transform_pair(z)
        half[k,:nr]=factor*z*p.speed*(m.nativeV[0]@dq);half[k,nr:]=factor*(m.nativeV[1]@dq)
        mult=(1 if k==go.COUNT//2 else 2)*kernel[k]*factor
        base_end+=(mult*q).real;corr_end+=(mult*dq).real;dEend+=(mult*dE).real
        if k%512==0:print('FEEDBACK',k,round(time.monotonic()-started,2),flush=True)
    half[-1]=half[-1].real;native=np.empty((3,65,2*nr),np.longdouble);coarse=np.empty((65,2*nr),np.longdouble)
    for lo in range(0,2*nr,128):
        hi=min(lo+128,2*nr);a=go.coefficients(half[:,lo:hi]);b=go.coefficients(half[::2,lo:hi])
        for j,n in enumerate(go.DEGREES):native[j,:,lo:hi]=L[:,:n]@a[:n]
        coarse[:,lo:hi]=L@b[:2048]
    initial_raw=float(abs(native[:,0]).max());native[:,0]=0;coarse[0]=0
    baseline=np.load(patch.OUT/'p4-2048.npz');w=baseline['weights'];masks=baseline['masks']
    def norm(values):
        v=values[...,:nr];f=values[...,nr:]
        return np.stack([np.sqrt(np.sum(w*v*v,axis=-1)),np.sqrt(np.sum(w*f*f,axis=-1)),
            *[np.sqrt(np.sum(w[mask]*v[...,mask]**2,axis=-1)/w[mask].sum()) for mask in masks]],axis=-1)
    base_native=np.c_[baseline['native_velocity'],baseline['native_scalar']]
    scale=norm(base_native).max(0); effect=norm(native[-1]).max(0)/scale
    time_error=norm(native[-1]-native[-2]).max(0)/scale;contour_error=norm(native[-1]-coarse).max(0)/scale
    endpoint_relative=float(abs(base_end-baseline['q']).max()/max(abs(baseline['q']).max(),1e-100))
    balance=float(abs(np.sum(-np.diff(dEend),dtype=np.longdouble))/max(abs(dEend).max(),1e-100))
    np.savez_compressed(OUT/'response.npz',times=go.TIMES,radius=m.original.native,weights=w,masks=masks,
        correction=native,coarse_correction=coarse,base_q=base_end,correction_q=corr_end,correction_heat_energy=dEend,
        total_native_velocity=baseline['native_velocity']+native[-1,:,:nr],total_native_scalar=baseline['native_scalar']+native[-1,:,nr:])
    passed=bool(max(time_error)<1e-6 and max(contour_error)<1e-6 and endpoint_relative<1e-10 and balance<2e-13 and p.loop_residual<1e-10)
    row=dict(classification='Counterexample candidate',actual_temperature_feedback_evolved=True,numerical_target_passed=passed,
        temperature_only_feedback_below_target=bool(passed and max(effect)<1e-4),
        comparisons={key:dict(correction_relative=float(effect[j]),time_absolute_to_baseline=float(time_error[j]),contour_absolute_to_baseline=float(contour_error[j])) for j,key in enumerate(go.task.FIELDS)},
        baseline_endpoint_replay=endpoint_relative,heat_balance=balance,linear_residual=p.error,feedback_residual=p.loop_residual,
        maximum_corrections=max(p.iterations,pilotrow['corrections']),sampled_second_to_first=max(p.second_to_first,pilotrow['second_to_first']),
        sampled_heat_feedback_relative=max(p.feedback_relative,pilotrow['sampled_heat_feedback_relative']),raw_initial_correction_max=initial_raw,
        seconds=time.monotonic()-started,total_compute_seconds=time.monotonic()-started+pilotrow['seconds'],memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        full_temperature_metric_coefficient_feedback=False,full_GR_photon_feedback_evolved=False,full_nonlinear_evolution=False,full_dynamic_charge_solved=False,
        scope='Temperature-only feedback on the fixed transport geometry/rates/coefficients. This new loop is solved at p4 with absolute correction tolerance relative to the accepted baseline, not a new spatial/physical certification or rigorous continuum error bound.')
    write(OUT/'result.json',row);signal.alarm(0);print('FEEDBACK RESULT',json.dumps(row),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run'])
    globals()[parser.parse_args().action]()
