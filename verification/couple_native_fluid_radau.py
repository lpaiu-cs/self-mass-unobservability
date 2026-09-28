"""Counterexample candidate: photon/E/H/B/S in the same actual Radau stages.

Reuse the native flux and photon owners. Finite local Jacobians propose Newton
steps; the original directional flux decides the accepted physical equation.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import gc,inspect,json,resource,shutil,sys,time
import numpy as np
import sympy as sp
from scipy import sparse
from scipy.sparse.linalg import LinearOperator,splu,gmres
import return_native_stage_collisions as history

prior=history.prior;previous=prior.previous;radau=prior.radau
OUT=Path('native-fluid-radau179-work');OLD=history.OUT
read,write,sha=prior.read,prior.write,prior.sha;LD=prior.LD;AMP=prior.AMP
A,B,C=radau.RK_A,radau.RK_B,radau.RK_C
CAPS=dict(prepare=15,check=80,pilot=480)
CAPS['check_retry']=50  # After20.862s rejected fixed-support check, within80s.
CAPS.update(branch_check=35,branch_pilot=420)
CAPS.update(branch64=240,branch128=420)


def selected_tangent():
    """Linear action on branches selected by one direction; never a new RHS."""
    decisions=[];position=0;recording=True
    def select(condition):
        nonlocal position
        if recording:decisions.append(np.array(condition,copy=True))
        answer=decisions[position];position+=1;return answer
    def minimum(a,b):return np.where(select(a<=b),a,b)
    def maximum(a,b):return np.where(select(a>=b),a,b)
    def minmod(a,b):return np.where(select(a*b>0),np.where(select(abs(a)<=abs(b)),a,b),0.)
    def clone(fn,changes,extra):
        source=inspect.getsource(fn)
        for old,new in changes:
            assert source.count(old)==1,(fn.__name__,old);source=source.replace(old,new)
        namespace=dict(fn.__globals__,**extra);exec(compile(source,__file__,'exec'),namespace);return namespace[fn.__name__]
    direction=clone(previous.engine.face.minmod_direction,[
        ('np.minimum(da,db)','minimum(da,db)'),('np.maximum(da,db)','maximum(da,db)'),
        ('(da*b>0)','select(da*b>0)'),('(db*a>0)','select(db*a>0)')],
        dict(minimum=minimum,maximum=maximum,minmod=minmod,select=select))
    reconstruction=clone(previous.engine.face.reconstruction,[('np.maximum(d[0],0)','maximum(d[0],0)')],
        dict(minmod_direction=direction,maximum=maximum))
    def extremum(values,derivatives,minimum):
        values=np.array(values);derivatives=np.array(derivatives)
        best=np.min(values,axis=0) if minimum else np.max(values,axis=0)
        masked=np.where(values==best,derivatives,np.inf if minimum else -np.inf)
        index=select(np.argmin(masked,axis=0) if minimum else np.argmax(masked,axis=0))
        return best,np.take_along_axis(derivatives,index[None],axis=0)[0]
    flux=FunctionType(previous.engine.flux_direction.__code__,dict(previous.engine.flux_direction.__globals__,extremum=extremum))
    deep=clone(previous.original.reuse.deep.deep_tangent,[('dmass>=0','select(dmass>=0)')],dict(select=select))
    # The owner tangent is generated; inspect its unexpanded source owner instead.
    source=previous.engine.source.replace('flux[0]>=0','select(flux[0]>=0)').replace('m.deep_tangent(k,z,field)','deep(m,k,z,field)')
    namespace=dict(previous.engine.tangent.__globals__,reconstruction=reconstruction,flux_direction=flux,deep=deep,select=select)
    exec(compile(source,__file__,'exec'),namespace);tangent=namespace['tangent']
    def apply(m,k,z,probe):return tangent(m,k,z,probe)
    def reset(capture):
        nonlocal position,recording
        if not recording:assert position==len(decisions),(position,len(decisions))
        position=0;recording=capture
    return apply,reset


def paths(sweep):return OUT/f'sweep-{sweep}/photons',OUT/f'sweep-{sweep}/material'


def prepare():
    assert not OUT.exists() and not read(OLD/'material-result.json')['passed'];OUT.mkdir()
    for s in [0,1]:
        for folder in paths(s):folder.mkdir(parents=True)
    for src in (OLD/'sweep-0').rglob('*.npz'):
        dst=OUT/src.relative_to(OLD);shutil.copyfile(src,dst);assert sha(src)==sha(dst)
    for name in ['normalization.json','photon-conservation-plan.json']:shutil.copyfile(OLD/name,OUT/name)
    files=[Path(__file__),Path(history.__file__),Path(prior.__file__),Path(previous.engine.__file__),
        Path(previous.engine.face.__file__),Path(radau.__file__),OLD/'material-result.json',OLD/'saved-state-check.json']
    files.extend(OLD/name for n in [64,128] for name in [f'material-{n}.npz',f'sweep-1/photons/pilot-{n}.npz'])
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='da36800b2',previous_turn='progress: actual stage collisions reduced the failed H transport below its original gate; baryon time and material-input closure remained rejected.',
        claim='Put B,S,Etilde,H and photons in one Radau equation, including actual native transport and photon momentum recoil, so B/S and mechanical energy/H are no longer lagged17-knot material inputs.',
        method='Keep the existing photon transport/collision, affine physical drive, native primitive/minmod/HLL/deep-gravity owners and Radau clocks. Four gas columns are Etilde/eu,H/nu,B/bu,S/su. Eref=Etilde+kappa*B. Native thermal rate is Eref_rate-kappa*B_rate. Inventory xi follows the same current B. No post-solve state overwrite.',
        proposal='A local finite Jacobian of the already-small native directional rate is a Newton proposal only. Individual deep/join columns include inventory nonlocality; five atmospheric colors resolve a five-cell rate stencil. Rebuild at solved directions if necessary, at most3 Newton solves. Original nonlinear native equation decides acceptance; no global linearity certificate.',
        residual='Require original1e-12 full equation and1e-13 physical E/H moment residuals. Also require1e-13 B/S physical residuals relative to RHS or stage-state size, with1e-290 floor; zero explicit B RHS must not erase the B equation. Keep original1e-8 balance,0.2percent constitutive/proposal tests and2percent time tests.',
        floor='Assert the active material support is fixed throughout the prefix. Inactive material RHS is projected to zero and its discarded native rate is recorded. A changing support stops this trial rather than silently deleting a stage state.',
        scope='Retained response with zero additional geometry. Original background, EOS and source failure records remain frozen. This is not full nonlinear GR or a final charge result.',
        reuse='Use178saved stage gas and free material as Newton initial directions only; no accepted unknown is fixed to them. Reuse old background and input arrays only to construct existing owners. Zero all lagged material sources before stages.',
        forecast='178pair measured147.924s; additional sparse native Jacobians and two material columns increase work. Check representative Jacobian timing first. One4/8macro pair,480s hard cap; assumed range150..480s, not guaranteed. No full period authorized.',
        budget=CAPS,CPU_threads=1,virtual_GiB=4,
        stop='Any control, original gate, Newton/GMRES limit, support change or cap. No automatic finer clock, method ladder, longer period or extra waveform sweep.',
        final_charge_conclusion='unadjudicated',bindings={str(p):sha(p) for p in files}))
    K,br,er=sp.symbols('K br er');assert sp.expand((er-K*br)+K*br-er)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='Fixed-background Eref=Etilde+kappa*B transforms native rates as Etilde_dot=Eref_dot-kappa*B_dot. It does not remove mass work or prove a full GR error bound.'))


def initialize():
    global Model
    previous.OUT=OUT;previous.paths=paths;previous.initialize();Parent=previous.run.c.Response
    feedback=previous.engine.face.feedback;two=feedback.mono.old.CoupledResponse
    class Coupled(Parent):
        def __init__(self,n):
            super().__init__(n)
            self.motion[:]=0;self.mechanical[:]=0;self.energy_offset[:]=0;self.xi[:]=0
            self.material.face_owner_error=0.
            self.material.deep_tangent=previous.original.reuse.deep.deep_tangent.__get__(self.material,type(self.material))
            self.kappa=(self.a.astype(LD)-self.model.m.a0)*self.model.cx*LD(previous.original.C)**2
            self.bu=self.eu/(self.model.cx*previous.original.C**2);self.su=self.eu.copy()
            self.units=np.stack([self.eu,self.nu,self.bu,self.su],axis=1)
            self.floor_discard=np.zeros((self.n,4),LD);self.floor_discard_history=[]
            self.stage_log=[];self.stage_t=[];self.stage_h=[];self.stage_states=[];self.stage_native=[];self.stage_collision=[];self.stage_discard=[]
            self.conserved_history=[];self.newton_iterations=[]
            self.guides=dict(np.load(OLD/f'sweep-1/photons/pilot-{n}.npz'))
            self.guide_material=dict(np.load(OLD/f'material-{n}.npz'));self.guide_C,self.guide_edges=history.dense(self.guides)
        def unpack(self,v):return v[:self.size].reshape(self.n,self.q,self.nf),v[self.size:].reshape(self.n,4)
        def conserved(self,g):
            z=np.array([g[:,2]*self.bu,g[:,3]*self.su,g[:,0]*self.eu,g[:,1]*self.nu],LD)
            z[2]+=self.kappa*z[0];return z
        def gas(self,p,b,e):
            eh=two.gas(self,p,b,e)
            recoil=-(np.sum(p*self.Eweight*self.mu[None,:,None],axis=(1,2))+e[2])/(self.a*self.su)
            return np.column_stack([eh,np.zeros(self.n),recoil])
        def inventory(self,g):return np.r_[self.material.model.mech.xi@np.r_[LD(0),-np.cumsum(g[:self.nb,2]*self.bu[:self.nb],dtype=LD)],np.zeros(self.n-self.nb)]
        def local(self,t):
            c=super().local(t);k=int(np.clip(np.searchsorted(self.t,t,side='right')-1,0,15));w=(t-self.t[k])/(self.t[k+1]-self.t[k])
            points=[(1-w,self.point(k)),(w,self.point(k+1))]
            for old,new in [('D','DBS'),('Db','DBSb'),('De','DBSe')]:
                c[new]=sum(v*r[old][...,:2]*np.stack([self.bu/r['units'][0],self.su/r['units'][1]],axis=-1).reshape(((1,self.n,2) if old=='De' else (self.n,1,1,2))) for v,r in points)
                c[new+'xi']=sum(v*r[old][...,3] for v,r in points)
            c['PBS']=sum(v*r['pressure_drive'][:,:2]*np.stack([self.bu/r['units'][0],self.su/r['units'][1]],axis=-1) for v,r in points)
            c['Pxi']=sum(v*r['pressure_drive'][:,3] for v,r in points);c['time']=t
            return c
        def collision(self,c,x,g,source=False):
            if g.shape[1]==2:
                assert not np.any(g),'Only the zero-gas lifted-source constructor may use two columns'
                g=np.column_stack([g,np.zeros_like(g)])
            eh=g[:,:2]
            p=-c['loss']*x+(c['S']@x.ravel()).reshape(x.shape)+np.einsum('nqfj,nj->nqf',c['B'],eh)
            b=-c['loss']*x+np.einsum('nqfj,nj->nqf',c['Bb'],eh)
            e=np.einsum('knqf,nqf->kn',c['esc'],x)+np.einsum('knj,nj->kn',c['Be'],eh)
            if 'DBS' in c:
                xi=self.inventory(g)
                p+=np.einsum('nqfj,nj->nqf',c['DBS'],g[:,2:])+c['DBSxi']*xi[:,None,None]
                b+=np.einsum('nqfj,nj->nqf',c['DBSb'],g[:,2:])+c['DBSbxi']*xi[:,None,None]
                e+=np.einsum('knj,nj->kn',c['DBSe'],g[:,2:])+c['DBSexi']*xi[None]
            if source:p+=c['q'];b+=c['qb'];e+=c['qe']
            return p,self.gas(p,b,e),e,b
        def native(self,t,g,probe=1.,details=False,tangent=None):
            z=self.conserved(g);k=int(np.clip(np.searchsorted(self.t,t,side='right')-1,0,15));w=(t-self.t[k])/(self.t[k+1]-self.t[k])
            F=np.zeros((4,self.n+1),LD);gravity=np.zeros(self.n,LD)
            for j,v in [(k,1-w),(k+1,w)]:
                if not v:continue
                f,r,_=(tangent or previous.engine.tangent)(self.material,j,z,probe);F+=v*f;gravity+=v*r
            raw=-np.diff(F,axis=1);raw[1]+=gravity;discard=np.zeros(4,LD)
            rate=raw
            normalized=np.column_stack([(rate[2]-self.kappa*rate[0])/self.eu,rate[3]/self.nu,rate[0]/self.bu,rate[1]/self.su])
            return (normalized,raw,discard,F,gravity) if details else normalized
        def jacobian(self,t,g):
            base=self.native(t,g);tangent,reset=selected_tangent()
            selected=self.native(t,g,tangent=tangent)
            assert np.array_equal(selected,base),'Selected branches must reproduce the original native direction'
            delta=np.maximum(np.max(abs(g),axis=0),LD('1e-100'))
            def action(q):
                reset(False);return self.native(t,q,tangent=tangent)
            rows=[];cols=[];values=[];join=self.nb+2
            for component in range(4):
                for j in range(join):
                    q=np.zeros_like(g);q[j,component]=delta[component];d=action(q)/delta[component]
                    ii=np.flatnonzero(np.any(d!=0,axis=1))
                    for r in range(4):
                        rows.extend(4*ii+r);cols.extend([4*j+component]*len(ii));values.extend(d[ii,r])
                for color in range(5):
                    jj=np.arange(join+color,self.n,5);q=np.zeros_like(g);q[jj,component]=delta[component]
                    d=action(q)
                    for offset in range(-2,3):
                        ii=jj+offset;ok=(ii>=0)&(ii<self.n);ri,cj=ii[ok],jj[ok]
                        for r in range(4):
                            rows.extend(4*ri+r);cols.extend(4*cj+component);values.extend(d[ri,r]/delta[component])
            J=sparse.coo_matrix((np.asarray(values,float),(rows,cols)),shape=(4*self.n,)*2).tocsr();J.eliminate_zeros()
            defect=(J@g.ravel()).reshape(self.n,4)-base
            relative=np.sum(abs(defect)*self.units,axis=0)/np.maximum(np.sum(abs(base)*self.units,axis=0),LD('1e-290'))
            assert max(relative)<1e-12,('Selected native branch reconstruction',relative.tolist())
            return J,base
        def inverse_pair(self,c,h,J):
            pair=SimpleNamespace(**vars(self));pair.pack=lambda x,g:np.r_[x.ravel(),g.ravel()]
            pair.unpack=lambda v:(v[:self.size].reshape(self.n,self.q,self.nf),v[self.size:].reshape(self.n,2))
            pair.gas=lambda p,b,e:two.gas(self,p,b,e)
            inv=two.inverse(pair,c,h);fluid=splu(sparse.eye(4*self.n,format='csc')-h*J)
            def apply(x,g):
                xx,eh=pair.unpack(inv(pair.pack(x,g[:,:2])));gg=g.copy();gg[:,:2]=eh
                return self.pack(xx,fluid.solve(np.asarray(gg.ravel(),float)).reshape(self.n,4))
            return apply
        def pressure(self,t,g):
            c=self.local(t)
            return np.einsum('nj,nj->n',c['pressure_map'],g[:,:2])+np.einsum('nj,nj->n',c['PBS'],g[:,2:])+c['Pxi']*self.inventory(g)
        def guide(self,t):
            ids=int(np.argmin(abs(self.guides['accepted_angular_times']-t)));assert abs(self.guides['accepted_angular_times'][ids]-t)<1e-18
            d=self.guide_material;j=int(np.clip(np.searchsorted(d['t'],t,side='left')-1,0,len(d['t'])-2));w=(t-d['t'][j])/(d['t'][j+1]-d['t'][j])
            z=(1-w)*d['lifted_history_scaled'][j]+w*d['lifted_history_scaled'][j+1]+self.guide_C(t)
            gas=self.guides['gas_stage_scaled'][ids]
            return np.column_stack([gas[:,0]/self.eu,gas[:,1]/self.nu,z[0]/self.bu,z[1]/self.su]).astype(LD)
    Model=Coupled
    source=(OUT/'sweep-1/expanded-corrected-run.py').read_text()
    def change(a,b):
        nonlocal source
        assert source.count(a)==1,(a,source.count(a));source=source.replace(a,b)
    change('g=np.zeros((self.n,2));','g=np.zeros((self.n,4));')
    change('transfer=np.zeros((self.n,2));','transfer=np.zeros((self.n,2),dtype=LD);')
    change("np.einsum('nj,nj->n',pmap,g)","self.pressure(t,g)")
    change('np.array([np.interp(t,self.t,v) for v in self.energy_offset.T])','self.kappa*g[:,2]*self.bu')
    change('np.array([np.interp(count*base_h,self.t,v) for v in self.energy_offset.T])','self.kappa*g[:,2]*self.bu')
    change('(1-gamma)*g1+gamma*g2-mechanical','(1-gamma)*g1[:,:2]+gamma*g2[:,:2]-mechanical')
    change('        physical_x=x+self.lift(t)[0]','        self.conserved_history.append(self.conserved(g)*AMPLITUDE)\n        physical_x=x+self.lift(t)[0]')
    change('            x,g=z,gz;moment=max(moment,e1,e2)',
        '            x,g=z,gz;moment=max(moment,e1,e2)\n            active=self.material.active(t+h);removed=g.copy();removed[active]=0.;g=g.copy();g[~active]=0.\n            self.floor_discard+=removed*self.units\n            ledger+=np.array([np.sum(removed[:,1]*self.nu),-np.sum(removed[:,0]*self.eu)],dtype=LD)')
    change('        self.conserved_history.append(self.conserved(g)*AMPLITUDE)',
        '        self.floor_discard_history.append(self.floor_discard.copy())\n        self.conserved_history.append(self.conserved(g)*AMPLITUDE)')
    change('energy_offset_reference=AMPLITUDE*self.energy_offset,energy_offset_t=self.t,',
        'energy_offset_reference=np.array(self.conserved_history)[:,0]*self.kappa,energy_offset_t=times,conserved_material_history=self.conserved_history,material_floor_discard_scaled=self.floor_discard,material_floor_discard_history_scaled=self.floor_discard_history,')
    change('additional_material_motion_evolved=False','additional_material_motion_evolved=True')
    ns=dict(Parent.run.__globals__,stages=stages,LD=LD)
    exec(compile(source,__file__,'exec'),ns);Coupled.run=ns['run'];(OUT/'expanded-joint-run.py').write_text(source)


def physical_norm(m,v):
    pairs=[m.unpack(row) for row in v.reshape(2,-1)];x=np.concatenate([p[0] for p in pairs]);g=np.concatenate([p[1] for p in pairs])
    return np.array([np.sum(abs(x)*np.tile(m.Nweight,(2,1,1)),dtype=LD)+np.sum(abs(g[:,1])*np.tile(m.nu,2),dtype=LD),
        np.sum(abs(x)*np.tile(m.Eweight,(2,1,1)),dtype=LD)+np.sum(abs(g[:,0])*np.tile(m.eu,2),dtype=LD),
        np.sum(abs(g[:,2])*np.tile(m.bu,2),dtype=LD),np.sum(abs(g[:,3])*np.tile(m.su,2),dtype=LD)])


def scales(m,rhs,sol):
    n=physical_norm(m,rhs);s=physical_norm(m,sol);n[:2]=np.maximum(n[:2],1.);n[2:]=np.maximum(n[2:],s[2:]);return np.maximum(n,LD('1e-290'))


def solve(m,op,P,rhs,guess):
    iterations=[];options=dict(M=P,rtol=1e-14,atol=0.,restart=20,maxiter=5,callback=iterations.append,callback_type='pr_norm')
    sol,info=gmres(op,np.asarray(rhs,float),x0=np.asarray(guess,float),**options);sol=sol.astype(LD)
    for k in range(4):
        residual=rhs-op.matvec(sol);relative=float(np.linalg.norm(residual)/max(np.linalg.norm(rhs),1e-290))
        moments=physical_norm(m,residual)/scales(m,rhs,sol)
        if relative<1e-14 and max(moments)<1e-13:
            radau.prior.owner.reuse.LINEAR.append(dict(initial_info=int(info),corrections=k,extended_residual=relative))
            m.max_residual=max(m.max_residual,relative);m.max_iterations=max(m.max_iterations,len(iterations));return sol
        assert k<3,('Four-moment linear residual',relative,moments.tolist())
        delta,_=gmres(op,np.asarray(residual,float),x0=None,**dict(options,rtol=1e-12));sol+=delta


def stages(m,t,h,x,g,lus):
    cs=[m.local(t+c*h) for c in C];ss=[m.source(t+c*h) for c in C];v=m.pack(x,g);dim=len(v)
    guides=[m.guide(t+c*h) for c in C];guess=np.array([m.pack(x,gg) for gg in guides]).ravel();audit=[]
    for newton in range(3):
        maps=[m.jacobian(t+c*h,gg) for c,gg in zip(C,guides)];Js=[r[0] for r in maps]
        affine=[r[1]-(J@gg.ravel()).reshape(m.n,4) for r,J,gg in zip(maps,Js,guides)]
        inv=[m.inverse_pair(c,h*A[j,j],Js[j]) for j,c in enumerate(cs)]
        def L(j,value):
            xx,gg=m.unpack(value);p,q,*_=m.collision(cs[j],xx,gg)
            stream=(m.A@xx.reshape(m.n*m.q,m.nf)).reshape(xx.shape)
            return m.pack(stream+p,q+(Js[j]@gg.ravel()).reshape(m.n,4))
        def mat(value):
            vv=value.reshape(2,dim);return (vv-h*(A@np.array([L(j,row) for j,row in enumerate(vv)]))).ravel()
        def pre(value):
            rows=[]
            for j,row in enumerate(value.reshape(2,dim)):
                xx,gg=m.unpack(row);xx=lus[j].solve(xx.reshape(m.n*m.q,m.nf)).reshape(xx.shape);rows.append(inv[j](xx,gg))
            return np.array(rows).ravel()
        src=np.array([m.pack(s[0]/(m.scale*AMP)+c['q'],m.gas(c['q'],c['qb'],c['qe'])+a) for s,c,a in zip(ss,cs,affine)])
        rhs=(np.tile(v,(2,1))+h*(A@src)).ravel();op=LinearOperator((2*dim,)*2,mat,dtype=float);P=LinearOperator(op.shape,pre,dtype=float)
        sol=solve(m,op,P,rhs,guess);pairs=[m.unpack(row) for row in sol.reshape(2,dim)];rates=[];details=[]
        for j,(xx,gg) in enumerate(pairs):
            p,q,e,b=m.collision(cs[j],xx,gg,True);native=m.native(t+C[j]*h,gg,details=True)
            rates.append(m.pack((m.A@xx.reshape(m.n*m.q,m.nf)).reshape(xx.shape)+p+ss[j][0]/(m.scale*AMP),q+native[0]))
            details.append((p,q,e,native))
        defect=(sol.reshape(2,dim)-v-h*(A@np.array(rates))).ravel();relative=float(np.linalg.norm(defect)/max(np.linalg.norm(rhs),1e-290))
        moments=physical_norm(m,defect)/scales(m,rhs,sol);audit.append(dict(relative=relative,moments=moments.astype(float).tolist()))
        if relative<1e-12 and max(moments)<1e-13:break
        if newton==2:
            np.savez_compressed(OUT/'rejected-joint-stage.npz',time=t,step=h,initial=v,solution=sol,
                guides=guides,defect=defect,physical_scales=scales(m,rhs,sol),native_rates=[d[3][0] for d in details])
            write(OUT/'rejected-joint-stage.json',dict(classification='Counterexample candidate',time=float(t),step=float(h),equations=audit))
            raise AssertionError(('True native joint Radau equation',audit))
        guess=sol;guides=[gg.copy() for _,gg in pairs]
    m.newton_iterations.append(audit);result=[];transport=[]
    for j,((xx,gg),(p,q,e,native)) in enumerate(zip(pairs,details)):
        now=t+C[j]*h;nr,raw,discard,F,gravity=native
        probes=[m.native(now,gg,v) for v in [.5,2.]]
        den=np.maximum(np.sum(abs(nr)*m.units,axis=0,dtype=LD),LD('1e-290'))
        probe=[(np.sum(abs(v-nr)*m.units,axis=0,dtype=LD)/den).astype(float).tolist() for v in probes]
        assert np.max(probe)<.002,('Native constitutive control',now,probe)
        native_balance=np.sum(raw,axis=1,dtype=LD)-(F[:,0]-F[:,-1]+np.r_[LD(0),np.sum(gravity,dtype=LD),LD(0),LD(0)])
        balance=(abs(native_balance)/np.maximum(np.sum(abs(raw),axis=1,dtype=LD),LD('1e-290'))).astype(float).tolist()
        assert max(balance)<1e-8,('Native shared-face balance',balance)
        m.stage_log.append(dict(time=float(now),probe=probe,native_balance=balance))
        m.stage_t.append(now);m.stage_h.append(h*B[j]);m.stage_states.append(m.conserved(gg));m.stage_native.append(nr*m.units);m.stage_collision.append(q*m.units);m.stage_discard.append(discard)
        result.append((xx,gg,p,q+nr,e,ss[j][1]/AMP,ss[j][2]));transport.append(nr[:,:2])
    return result,sum(b*r for b,r in zip(B,transport))


def check(retry=False,branch=False):
    if retry:
        failure=read(OUT/'check-receipt.json');assert failure['error']=="AssertionError('Changing material support')"
        assert failure['seconds']+CAPS['check_retry']<80
        write(OUT/'floor-repair-plan.json',dict(classification='Counterexample candidate',failure=failure,
            change='The background activity mask changes inside the prefix. Withdraw fixed-support assumption. Restore the existing native material floor projection after each accepted Radau substep, applying its actual endpoint mask to all four gas variables. Record removed conserved quantities, subtract their Etilde from the energy ledger and add their H to the photon-minus-H ledger. No stage solution is replaced by another trajectory.',
            physical_equation='Native flux and collisions evolve all cells in both stages, as in the original free material RHS. The floor is a separate explicit finite map at accepted boundaries. Test the true pre-projection equations and post-projection balances, and store both stage and projected endpoint states.',
            original_check_seconds=80,remaining_check_seconds=CAPS['check_retry'],full_horizon_authorized=False,
            bindings={str(p):sha(p) for p in [Path(__file__),OUT/'initial-producer.py',OUT/'fixed-support-plan.json',OUT/'check-receipt.json']}))
    initialize();m=Model(64);now=m.guides['accepted_angular_times'][3];g=m.guide(now);started=time.monotonic();J,base=m.jacobian(now,g);seconds=time.monotonic()-started
    phase=np.arange(m.n)[:,None]+np.arange(4)[None];direction=g*np.cos(phase*.7);h=LD(2)**-12
    exact=(m.native(now,g+h*direction)-base)/h;proposal=(J@direction.ravel()).reshape(m.n,4)
    den=np.maximum(np.sum(abs(exact)*m.units,axis=0,dtype=LD),LD('1e-290'));errors=(np.sum(abs(exact-proposal)*m.units,axis=0,dtype=LD)/den).astype(float).tolist()
    assert max(errors)<.002,('Local native Jacobian',errors)
    c=m.local(now);zero=np.zeros_like(m.I[0]);p,q,e,b=m.collision(c,zero,g)
    # Direct owner D mapping, including B-induced inventory, must give the same photon source.
    k=int(np.clip(np.searchsorted(m.t,now,side='right')-1,0,15));w=(now-m.t[k])/(m.t[k+1]-m.t[k]);xi=m.inventory(g)
    target=np.zeros_like(p)
    for weight,row in [(1-w,m.point(k)),(w,m.point(k+1))]:
        drive=np.stack([g[:,2]*m.bu/row['units'][0],g[:,3]*m.su/row['units'][1],np.zeros(m.n),xi],axis=-1)
        target+=weight*(np.einsum('nqfj,nj->nqf',row['D'],drive)+np.einsum('nqfj,nj->nqf',row['B'],g[:,:2]))
    mapping=float(np.sum(abs(target-p)*m.Eweight,dtype=LD)/max(np.sum(abs(target)*m.Eweight,dtype=LD),1.));assert mapping<1e-12
    forecast=read(OLD/'capture-receipt.json')['seconds']+72*seconds+30
    result=dict(classification='Counterexample candidate',passed=True,jacobian_seconds=seconds,jacobian_nnz=J.nnz,
        jacobian_directional_relative=errors,photon_current_material_mapping_relative=mapping,
        estimated_pilot_seconds=forecast,eligible=forecast<CAPS['branch_pilot' if branch else 'pilot'],
        scope='One actual saved direction; Jacobian is only a solver proposal. The trial must satisfy true native equations at every accepted stage.')
    write(OUT/('branch-check-result.json' if branch else 'check-result.json'),result);print(json.dumps(result),flush=True)


def branch_check():
    spent=sum(read(OUT/f'{name}-receipt.json')['seconds'] for name in ['check','check_retry'])
    assert spent+CAPS['branch_check']<CAPS['check']
    assert read(OUT/'pilot-receipt.json')['seconds']+CAPS['branch_pilot']<CAPS['pilot']
    write(OUT/'branch-plan.json',dict(classification='Conjectural',
        failure='The first joint stage stopped after three solves: S equation6.22e-10 exceeds1e-13 although the vector residual was2.31e-16. No stage was accepted.',
        root='The Jacobian subtracts nearly equal directional rates through a primitive map with float64 outputs; that can corrupt small momentum columns. Replace that subtraction by direct column actions with the native minmod, density, wave-speed and donor branches selected at the current guide. Do not alter the original physical RHS or its accepted precision.',
        controls='Selected action at the guide must exactly reproduce native RHS; J times guide must match all four native rate channels below1e-12. Compare a separate local direction below0.2percent. True native Radau equations still decide acceptance with original1e-12/1e-13 gates and max3Newton solves.',
        budget=dict(original_check=80,check_spent=spent,branch_check=CAPS['branch_check'],original_pilot=480,pilot_spent=read(OUT/'pilot-receipt.json')['seconds'],branch_pilot=CAPS['branch_pilot']),
        forecast='Original147.924s photon pair plus72measured Jacobians plus30s reserve;72 assumes two Newton solves per each of18physical steps, not a guaranteed bound. Hard wall cap420s remains. No new clock or full period.',
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'finite-jacobian-producer.py',OUT/'pilot-receipt.json',OUT/'floor-repair-plan.json']}))
    check(branch=True)


def branch_pilot():
    assert read(OUT/'branch-check-result.json')['eligible'];pilot()


def pilot(selected=(64,128)):
    assert read(OUT/'check-result.json')['eligible'];initialize();rows=[];ph=[];mh=[]
    for n in selected:
        m=Model(n);row=m.run(n,f'pilot-{n}',n//16);file=paths(1)[0]/f'pilot-{n}.npz';p=dict(np.load(file))
        p.update(joint_stage_times=np.array(m.stage_t),joint_stage_weights=np.array(m.stage_h),joint_stage_conserved_scaled=np.array(m.stage_states),
            joint_native_rates_scaled=np.array(m.stage_native),joint_collision_rates_scaled=np.array(m.stage_collision),joint_discard_rates_scaled=np.array(m.stage_discard))
        assert np.array_equal(p['joint_stage_times'],p['accepted_angular_times'])
        assert np.array_equal(p['joint_stage_weights'],p['accepted_angular_quadrature_weights'])
        weighted=np.sum(p['joint_stage_weights'][:,None,None].astype(LD)*(p['joint_native_rates_scaled']+p['joint_collision_rates_scaled']),axis=0,dtype=LD)
        gas=p['delta_material']/AMP;expected=gas*m.units
        residual=(np.sum(abs(expected+p['material_floor_discard_scaled']-weighted),axis=0,dtype=LD)/np.maximum(np.sum(abs(expected),axis=0,dtype=LD),LD('1e-290'))).astype(float).tolist()
        assert max(residual)<1e-8,('Same-solution gas ledger',residual)
        row.update(joint_material_ledger_relative=residual,actual_B_S_E_H_joint_unknowns=True,lagged_material_input_used=False,
            native_true_equation=max(v[-1]['relative'] for v in m.newton_iterations),native_true_moments=max(max(v[-1]['moments']) for v in m.newton_iterations),
            maximum_Newton_solves=max(len(v) for v in m.newton_iterations))
        np.savez_compressed(file,**p);write(file.with_suffix('.json'),row);write(OUT/f'stages-{n}.json',dict(classification='Counterexample candidate',equations=m.newton_iterations,controls=m.stage_log))
        rows.append(row);ph.append(p['moments'][:,[0,1,2,3,5,6]]);mh.append(p['conserved_material_history']);assert row['passed'],row
        del m,p;gc.collect()
    if selected!=(64,128):return rows
    return compare_paths()


def compare_paths():
    rows=[];ph=[];mh=[]
    for n in [64,128]:
        file=paths(1)[0]/f'pilot-{n}.npz';row=read(file.with_suffix('.json'));assert row['passed']
        with np.load(file) as p:ph.append(p['moments'][:,[0,1,2,3,5,6]]);mh.append(p['conserved_material_history'])
        rows.append(row)
    photon=previous.run.c.relative(*ph);material=previous.run.c.relative(*mh)
    result=dict(classification='Counterexample candidate',passed=max(photon+material)<.02,rows=rows,photon_time=photon,material_time=material,
        full_horizon_completed=False,GR_charge_readout_executed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'pilot-result.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


def branch64():
    assert read(OUT/'branch-check-result.json')['passed']
    write(OUT/'sequential-branch-plan.json',dict(classification='Conjectural',
        admission='Paired forecast425.685s exceeded the conservative420s allocation. Do not launch that pair blindly. Run only the original64prefix first, max240s; save its actual joint stages for reuse. The original480s budget includes the prior54.097s failure. Admit the original128prefix only if its measured forecast fits the remaining original budget.',
        remaining_original_pilot_seconds=480-read(OUT/'pilot-receipt.json')['seconds'],
        fine_forecast='Coarse stepping time times2 for twice the actual steps, plus measured coarse operator construction and20s I/O reserve. Late step/iteration cost is still an assumption. Strict remaining wall cap; no extra clock, period or rerun of accepted coarse prefix.',
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'finite-jacobian-producer.py',OUT/'pilot-receipt.json',OUT/'branch-check-result.json',OUT/'branch-plan.json']}))
    pilot((64,))


def branch128():
    coarse=read(paths(1)[0]/'pilot-64.json');spent=sum(read(OUT/f'{name}-receipt.json')['seconds'] for name in ['pilot','branch64'])
    forecast=2*coarse['stepping_seconds']+coarse['operator_point_seconds']+20
    write(OUT/'fine-admission.json',dict(classification='Counterexample candidate',spent=spent,remaining=480-spent,
        measured_coarse_seconds=coarse['seconds'],fine_forecast_seconds=forecast,eligible=forecast<480-spent,
        full_horizon_authorized=False,final_charge_conclusion='unadjudicated'))
    assert forecast<480-spent,('Original pilot budget admission',forecast,480-spent)
    pilot((128,));compare_paths()


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    cap=CAPS[action]
    if action=='branch128':cap=min(cap,int(480-sum(read(OUT/f'{name}-receipt.json')['seconds'] for name in ['pilot','branch64'])))
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));previous.original.inf.incident.native.deadline(cap)
    started=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():
                target=OUT/'initial-producer.py' if Path(p)==Path(__file__) and (OUT/'initial-producer.py').exists() else p
                assert sha(target)==h,p
            if action=='pilot':
                for p,h in read(OUT/'floor-repair-plan.json')['bindings'].items():assert sha(p)==h,p
            if action=='branch_pilot':
                for p,h in read(OUT/'branch-plan.json')['bindings'].items():assert sha(p)==h,p
            if action=='branch128':
                for p,h in read(OUT/'sequential-branch-plan.json')['bindings'].items():assert sha(p)==h,p
        if action=='check_retry':check(True)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-started,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
