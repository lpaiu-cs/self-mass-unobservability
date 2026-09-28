"""Actual orbital scalar drive through the coupled conductive GR response.

Use physical GR coordinates q and integrated active face energy E. The sparse
joint solve closes the temperature/current loop without a low-frequency
fixed-point assumption. Unknown outer transport remains an omitted component.
"""
from pathlib import Path
import argparse
import json
import signal
import resource
import time
import numpy as np
from scipy.sparse import bmat,diags,csc_matrix
from scipy.sparse.linalg import splu
import def_gr_temperature_feedback as thermal
import def_orbital_charge_fem as orbit

go=thermal.go;write=go.write
OUT=orbit.OUT.parent/'def-orbital-conductive-feedback'


class Problem(thermal.Problem):
    def __init__(self,degree=4,coarse=False):
        thermal.patch.install();super().__init__(degree);m=self.model;heat=m.heat
        self.active=heat.face_ids;self.H=m.H[:,self.active].astype(np.clongdouble)
        self.load_active=self.load[:,self.active]
        self.Gq=self.Gq.astype(np.clongdouble)
        self.GEraw=(self.GE+self.Gq@m.H)[:,self.active].astype(np.clongdouble)
        # A constant scalar is not in the fixed-zero outer trial space.
        a=m.data[:,:4].reshape(-1,2,2).astype(np.longdouble)
        b=m.data[:,4:8].reshape(-1,2,2).astype(np.longdouble)
        c=m.data[:,8:12].reshape(-1,2,2).astype(np.longdouble)
        w=m.data[:,12:16].reshape(-1,2,2).astype(np.longdouble);weights=m.weights.astype(np.longdouble)
        V=[v.astype(np.longdouble) for v in m.V];cov=[v.astype(np.longdouble) for v in m.cov];lift=-a[:,:,1]
        self.kl=sum(cov[i].T@(weights*sum(b[:,i,j]*lift[:,j] for j in range(2)))+V[i].T@(weights*c[:,i,1]) for i in range(2))
        self.ml=sum(V[i].T@(weights*w[:,i,1]) for i in range(2))
        self.kll=np.sum(weights*(np.einsum('ni,nij,nj->n',lift,b,lift)+c[:,1,1]),dtype=np.longdouble)
        self.mll=np.sum(weights*w[:,1,1],dtype=np.longdouble)
        source=go.task.fem.source_points(heat,m.points)
        gs=m.data[:,16:22].reshape(-1,2,3);hs=m.data[:,22:28].reshape(-1,2,3)
        g=[sum(diags(gs[:,i,j])@source[j] for j in range(3)) for i in range(2)]
        h=sum(diags(hs[:,1,j])@source[j] for j in range(3))
        fL=sum(g[j].T@(weights*sum(lift[:,i]*b[:,i,j] for i in range(2))) for j in range(2))-h.T@weights
        self.fL=np.asarray(fL).ravel()[self.active].astype(np.clongdouble)
        # The heat momentum lift has no scalar component and W is diagonal.
        assert float(abs(self.ml@self.H).max())==0
        point=self.point;r=point['r'];alpha=-4*point['phi']
        TL=-point['adiabatic_T_rho']*(r*point['v']+3*alpha)
        bank=np.load(go.task.BANK/'fine-bank.npz');faces=bank['faces'];n=len(r)
        self.gL=(self.theta*TL)[n-1-faces]-(self.theta*TL)[n-faces]
        self.TL=TL
        if coarse:
            other=np.load(go.task.BANK/'coarse-bank.npz');assert np.array_equal(faces,other['faces'])
            pref=self.conductance.sum(1)/bank['mode_K_SI'].sum(1)
            self.conductance=pref[:,None]*other['mode_K_SI']
            d=heat.d;AN=np.sqrt(d['A'][faces-1]*d['A'][faces]*d['N'][faces-1]*d['N'][faces])
            self.lam=other['poles_proper_s']*AN[:,None]*heat.geometry.tc
        # Energy scaling by the local thermal response, followed by sparse
        # row/column equilibration; ordering follows actual radial support.
        self.energy_scale=1/np.maximum(abs(self.GEraw.diagonal()),np.longdouble('1e-100'))
        positions=np.empty(m.size)
        for field in range(2):
            valid=m.indices[:,field]>=0;positions[m.indices[valid,field]]=m.grid[valid]
        self.permutation=np.argsort(np.r_[positions,heat.edges[self.active]],kind='stable')
        bg=m.bg.sample(np.array([2.]));self.mu=float(bg['m'][0]/2);self.flux=float(bg['v'][0]*2)
        self.F=float(bg['N'][0]*np.sqrt(1-2*self.mu));self.R=heat.geometry.R
        old=json.loads((orbit.OUT/'pilot.json').read_text());self.mass=old['ADM_geom_m']
        self.omega=old['row']['omega_R_over_c']
        from scipy.integrate import quad
        self.delay=quad(lambda x:float(heat.geometry.metric(np.array([x]))[1][0]/heat.geometry.metric(np.array([x]))[0][0]),1,2,epsabs=1e-12)[0]

    def solve(self,n,label):
        start=time.monotonic();m=self.model;z=np.clongdouble(-1j*n*self.omega);w=float(-z.imag)
        G=(self.K+complex(z*z)*self.M).tocsc();Gx=self.Kx+z*z*self.Mx
        rhs=-(self.kl+z*z*self.ml).astype(np.clongdouble)
        sc=np.sqrt(abs(G.diagonal()));D=diags(1/sc);lu=splu((D@G@D).astype(complex).tocsc(),permc_spec='NATURAL')
        qa=(lu.solve(np.asarray(rhs/sc,complex))/sc).astype(np.clongdouble)
        for _ in range(3):qa+=(lu.solve(np.asarray((rhs-self.Kx@qa-z*z*(self.Mx@qa))/sc,complex))/sc).astype(np.clongdouble)
        qa_error=float(np.max(abs(rhs-self.Kx@qa-z*z*(self.Mx@qa))/(abs(rhs)+self.absK@abs(qa)+abs(z*z)*(self.absM@abs(qa))+1e-100)))
        gain=np.sum(self.conductance*m.heat.geometry.tc*self.lam/(z*(z+self.lam)),axis=1)
        B=self.load_active-self.Kx@self.H-z*z*(self.Mx@self.H)
        C=diags(1/gain)-self.GEraw
        block=bmat([[Gx,-B@diags(self.energy_scale)],[-self.Gq,C@diags(self.energy_scale)]],format='csc')
        forcing=self.Gq@qa+self.gL
        bvec=np.r_[np.zeros(m.size,np.clongdouble),forcing]
        perm=self.permutation;A=block[perm,:][:,perm].astype(complex).tocsc();bp=bvec[perm]
        row=np.asarray(abs(A).max(axis=1).toarray()).ravel();Ar=diags(1/row)@A
        col=np.asarray(abs(Ar).max(axis=0).toarray()).ravel();scaled=(Ar@diags(1/col)).tocsc()
        factor=splu(scaled,permc_spec='NATURAL')
        def invert(b):
            yp=factor.solve(np.asarray(b[perm]/row,complex))/col
            y=np.empty_like(yp);y[perm]=yp
            return y.astype(np.clongdouble)
        answer=invert(bvec)
        for _ in range(4):answer+=invert(bvec-block@answer)
        residual=bvec-block@answer
        error=float(np.max(abs(residual)/(abs(bvec)+abs(block)@abs(answer)+1e-100)))
        dq=answer[:m.size];E=answer[m.size:]*self.energy_scale
        heat_drive=self.Gq@(qa+dq)+self.gL+self.GEraw@E
        heat_defect=E-gain*heat_drive
        heat_error=float(abs(heat_defect).max()/max(abs(E).max(),1e-100))
        # Full physical-coordinate and heat equations, independent of the
        # assembled correction block and its row/column scaling.
        gr_defect=self.Kx@dq+z*z*(self.Mx@dq)-B@E
        gr_error=float(np.max(abs(gr_defect)/(self.absK@abs(dq)+abs(z*z)*(self.absM@abs(dq))+abs(B)@abs(E)+1e-100)))
        assert max(error,qa_error,gr_error)<1e-9 and heat_error<1e-9,(error,qa_error,gr_error,heat_error)
        Ra=-rhs@qa+self.kll+z*z*self.mll
        dR=-rhs@dq-self.fL@E
        Za=Ra/(2*self.F);dZ=dR/(2*self.F)
        wave=orbit.exterior.outgoing(self.mu,self.flux,2*w);h=wave['h'];Zo=wave['impedance']
        drive=np.exp(-1j*w*(1+self.delay))/(self.F*h);Da=Za-Zo
        aa=drive/Da;af=drive/(Da+dZ);da=-drive*dZ/(Da*(Da+dZ))
        tail=2*da/h*np.exp(-1j*w*self.delay);charge=-self.R/self.mass*tail
        fullE=np.zeros(len(m.heat.edges),np.clongdouble);fullE[self.active]=af*E
        temp=af*(self.Tq@(qa+dq)+self.TL)+self.TE@fullE
        contrast_temp=af*(self.Tq@dq)+da*(self.Tq@qa+self.TL)+self.TE@fullE
        balance=float(abs(np.sum(-np.diff(fullE),dtype=np.clongdouble))/max(abs(fullE).max(),1e-100))
        # Real dissipated quadratic form of the same positive thermal poles.
        flux=z*E/m.heat.geometry.tc
        dissipation=float(np.real(np.vdot(heat_drive,flux)))
        assert dissipation>=0 and balance<2e-13
        old=json.loads((orbit.OUT/'result.json').read_text());incident=old['rows'][n]['drive_amplitude']
        pair=lambda v:[float(v.real),float(v.imag)]
        np.savez_compressed(OUT/f'{label}-{n}.npz',qa=qa,dq=dq,E=E,temperature=temp,temperature_correction=contrast_temp,
            active_faces=self.active,native_radius=m.original.native,grid=m.grid,indices=m.indices)
        result=dict(classification='Counterexample candidate',harmonic=n,degree=m.degree,seconds=time.monotonic()-start,
            scalar_boundary_value=pair(af),adiabatic_boundary_value=pair(aa),adiabatic_impedance=pair(Za),thermal_impedance_difference=pair(dZ),
            outgoing_conduction_contribution=pair(tail),radiative_charge_gain=pair(charge),
            actual_drive_amplitude=incident,delta_alpha_radiative_over_phi0=float(abs(charge)*incident/.001),
            maximum_actual_Eulerian_delta_lnT=float(abs(temp).max()*incident),
            maximum_actual_thermal_temperature_correction=float(abs(contrast_temp).max()*incident),
            linear_residual=error,adiabatic_residual=qa_error,physical_GR_residual=gr_error,heat_law_residual=heat_error,
            heat_balance=balance,positive_pole_dissipation=dissipation,
            active_faces=len(self.active),dofs=block.shape[0],matrix_nonzeros=block.nnz,
            factor_nonzeros=factor.L.nnz+factor.U.nnz,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
            full_radial_photons=False,full_nonlinear_or_thermal_background=False,full_goal_complete=False)
        write(OUT/f'{label}-{n}.json',result);print('ORBITAL HEAT',label,n,json.dumps(result),flush=True)
        return result


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Apply the actual declared orbital scalar drive to the full coupled conductive GR temperature response on all4012 computed internal faces, and extract the same-incoming outgoing thermal contribution.',
        model='Physical q and integrated active face energy E; (K+z^2 M)q=(F-z^2 M H)E+external drive; [diag(1/gain)-GEraw]E=Gq*q+gL. Original positive microscopic poles, native EOS thermal tangent, heat energy debit, mass constraint and heat momentum retained. Coefficients/geometry stay frozen.',
        scope='Driven harmonic perturbation relative to the same affine frozen background. This does not prove that the actual nonstationary background is fixed over an orbit. Unknown outer faces and photons are omitted contributions, never a certified insulating boundary.',
        cases='p4 and p2 harmonics1,2,3 plus original coarse-bank coefficient contrast. No zero-frequency thermal inverse or extra carriers.',
        gates=dict(linear_residual=1e-9,heat_law_residual=1e-9,heat_balance=2e-13,spatial=.02,coefficient=.02),
        budget=dict(pilot_seconds=120,total_compute_seconds=600,CPU_threads=1,memory_GB=4,new_EOS_calls=0,new_time_steps=0),
        decision='Reuse pilot. Stop on failed solve or spatial/coefficient gate, without automatic enlargement. Use solved heat correction, not fixed-point convergence extrapolated from short-time contours.',
        bindings={str(p):go.task.digest(p) for p in [Path(__file__),Path(thermal.__file__),Path(orbit.__file__),
            go.task.BANK/'fine-bank.npz',go.task.BANK/'coarse-bank.npz',orbit.OUT/'result.json',orbit.OUT/'pilot.json']}))
    write(OUT/'symbolic.json',thermal.control());signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)))
    started=time.monotonic();p=Problem();setup=time.monotonic()-started;row=p.solve(1,'p4')
    forecast=1.5*(3*setup+9*row['seconds']+20)
    write(OUT/'pilot.json',dict(classification='Counterexample candidate',setup_seconds=setup,row=row,forecast_seconds=forecast,
        seconds=time.monotonic()-started,assumption='Three p4-sized assemblies and nine joint solves scaled from measured pilot plus20s and50 percent margin; p2 and coarse-bank unmeasured.'))
    print('FORECAST',forecast,flush=True)


def run():
    plan=json.loads((OUT/'plan.json').read_text());pilot=json.loads((OUT/'pilot.json').read_text())
    for p,h in plan['bindings'].items():assert go.task.digest(Path(p))==h,p
    assert not (OUT/'result.json').exists() and pilot['forecast_seconds']<600
    signal.alarm(int(600-pilot['seconds']));resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));start=time.monotonic();cases={}
    for label,degree,coarse in [('p4',4,False),('p2',2,False),('coefficient',4,True)]:
        p=Problem(degree,coarse);cases[label]=[pilot['row'] if label=='p4' and n==1 else p.solve(n,label) for n in [1,2,3]];del p
        if label=='p2':
            spatial=[abs(complex(*a['radiative_charge_gain'])-complex(*b['radiative_charge_gain']))/abs(complex(*b['radiative_charge_gain'])) for a,b in zip(cases['p2'],cases['p4'])]
            if max(spatial)>.02:break
    comparisons=[]
    for i in range(3):
        ref=complex(*cases['p4'][i]['radiative_charge_gain']);item=dict(harmonic=i+1)
        for label in ['p2','coefficient']:
            if label in cases:item[label]=abs(complex(*cases[label][i]['radiative_charge_gain'])-ref)/abs(ref)
        comparisons.append(item)
    passed=all(r.get('p2',1)<.02 and r.get('coefficient',1)<.02 for r in comparisons)
    result=dict(classification='Counterexample candidate',passed=passed,cases=cases,comparisons=comparisons,
        seconds=time.monotonic()-start,total_compute_seconds=time.monotonic()-start+pilot['seconds'],
        orbital_conductive_GR_feedback_solved=True,full_radial_photons=False,full_nonlinear_or_thermal_background=False,
        static_thermal_comparator_solved=False,observational_nuisance_applied=False,full_goal_complete=False)
    write(OUT/'result.json',result);print('RESULT',json.dumps({k:v for k,v in result.items() if k!='cases'}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run']);globals()[parser.parse_args().action]()
