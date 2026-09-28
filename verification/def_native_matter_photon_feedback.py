"""Return conserved material motion to actual photon/energy/H evolution.

Counterexample candidate: a waveform-relaxation sweep. Baryon, momentum and
noncollisional material transport are supplied by the completed material sweep;
photons and total conserved energy/H are simultaneous unknowns. Iterate matter
and GR afterwards; this single return sweep is not full closure.
"""
from pathlib import Path
from types import FunctionType
import inspect,textwrap,json,signal,sys,time
import numpy as np
import sympy as sp
from scipy import sparse
import def_native_monolithic_response as mono
import def_native_material_branch_response as matter
import verify_native_material_response as primitive_owner

old=mono.old;flow=old.flow;C=old.C;write=old.write;sha=old.sha;AMP=old.AMPLITUDE
OUT=matter.OUT.parent.parent/'def-native-matter-photon-feedback'
replace=matter.base.replace

# Reuse the audited direct conservative-to-primitive variation. Omit its
# readout-only finite probes; local inventory displacement is a separate input.
primitive_source=inspect.getsource(primitive_owner.pressure).split('    # Verify pressure')[0]
primitive_source=replace(primitive_source,'def pressure(m,k,z,field):','def primitive(m,k,z,field,bank):')
primitive_source=replace(primitive_source,"    bank=dict(np.load(run.base.photons.old.OUT/f'bank-{m.reference}/point-{k}.npz'))\n",'')
primitive_source=replace(primitive_source,'xi=model.mech.xi@dh','xi=np.zeros(nb)')
# Use Etilde=Eref-(a_ref-a_surface)*cx*c^2*B as the energy coordinate.
# This is an invertible conservative-variable change, not removal of mass work.
primitive_source=replace(primitive_source,"du=np.asarray((z[2,:nb].astype(LD)-(m.a[:nb].astype(LD)-model.m.a0)*model.cx*LD(C)**2*z[0,:nb])/(m.a[:nb]*q[0,:nb]),float)","du=np.asarray(z[2,:nb].astype(LD)/(m.a[:nb]*q[0,:nb]),float)")
primitive_source=replace(primitive_source,"RE=(z[2,ids].astype(LD)+LD(m.rest)*z[0,ids])/(m.a[ids]*m.V[ids])-E*s[ids]-(W2*Hr-pr)*DD-(W2*Hy-py)*Y", "W=np.sqrt(W2);wminus1=W2*v*v/(W+1);w2minus1=W2*v*v\n    K0=rho*model.cx*LD(C)**2*W*wminus1+rho*uu*W2+pp*w2minus1\n    Kr=rho*model.cx*LD(C)**2*W*wminus1+rho*(uu+bank['dr_u'][ids])*W2+pr*w2minus1\n    Kt=rho*bank['dt_u'][ids]*W2+pt*w2minus1;Ky=rho*bank['dy_u'][ids]*W2+py*w2minus1\n    Kv=rho*model.cx*LD(C)**2*v*W**3*(2*W-1)+2*(rho*uu+pp)*W2*W2*v\n    RE=z[2,ids].astype(LD)/(m.a[ids]*m.V[ids])-K0*s[ids]-Kr*DD-Ky*Y")
primitive_source=replace(primitive_source,"ET=W2*Ht-pt;EV=2*H*W2*W2*v-(W2*Hr-pr)*W2*v", "ET=Kt;EV=Kv-Kr*W2*v")
primitive_source+='    return dict(dr=dr,dv=dv,dt=dt,dy=dy,dp=dp)\n'
namespace=dict(vars(primitive_owner));exec(compile(primitive_source,__file__,'exec'),namespace);primitive=namespace['primitive']


class Response(mono.Response):
    def __init__(self,reference):
        super().__init__(reference);self.material=matter.Material(reference)
        z=np.load(matter.OUT/f'material-source-128-reference-{reference}.npz')
        self.motion=np.stack([z['delta_baryon_g'],z['delta_radial_momentum_c_erg'],z['delta_reference_material_energy_erg'],z['delta_neutral_number']],axis=1)/AMP
        assert self.model.cx==self.model.flow.eos.cx
        self.energy_offset=(self.a-self.model.m.a0)*self.model.cx*C*C*self.motion[:,0]
        self.motion[:,2]-=self.energy_offset
        self.mechanical=self.motion[:,[2,3]]-self.material.transfer[:,[2,3]]
        self.xi=np.array([self.material.model.mech.xi@np.r_[0.,-np.cumsum(v[0,:self.nb])] for v in self.motion])
        self.velocity_jet_error=0.;self.mapping_error=0.;self.map_checks=[]

    def point(self,k):
        if k in self.cache:return self.cache[k]
        start=time.monotonic();self.restore(k);c=dict(np.load(old.OUT/f'bank-{self.reference}/point-{k}.npz'))
        S,esc=self.scattering_matrix(c);X=self.I[k]/self.scale;scat=(S@X.ravel()).reshape(X.shape);escaped=np.einsum('knqf,nqf->kn',esc,X)
        em=c['emit']/self.scale;bound=em-c['loss']*X;coll=bound+scat
        derivative={};dbound={};desc={}
        for p in ['dr','dt','dy']:
            ratio=np.divide(c[p+'_sc'],c['sc'],out=np.zeros(self.n),where=c['sc']>0)
            dbound[p]=c[p+'_emit']/self.scale-c[p+'_loss']*X
            derivative[p]=dbound[p]+ratio[:,None,None]*scat;desc[p]=ratio[None]*escaped
        # A one-sided Richardson derivative stays on the nonzero velocity
        # remap branch. At beta=0, local energy/H do not induce delta beta.
        beta0=self.beta.copy();bulk0=self.bulk_beta.copy();sgn=np.where(c['beta']<0,-1.,1.);vjets=[]
        for h in [2e-7,1e-7]:
            eps=sgn*h;self.beta=beta0+eps[self.nb:];self.bulk_beta=bulk0+eps[:self.nb]
            cc=self.coefficients();ss,ee=self.scattering_matrix(cc)
            db=(cc['emit']/self.scale-cc['loss']*X-bound)/eps[:,None,None]
            dv=db+((ss-S)@X.ravel()).reshape(X.shape)/eps[:,None,None]
            de=np.einsum('knqf,nqf->kn',ee-esc,X)/eps[None]
            vjets.append((dv,db,de))
        self.beta=beta0;self.bulk_beta=bulk0
        derivative['dv'],dbound['dv'],desc['dv']=[2*b-a for a,b in zip(*vjets)]
        err=float(np.sum(abs(vjets[1][0]-vjets[0][0])*self.Eweight)/max(np.sum(abs(derivative['dv'])*self.Eweight),1.));self.velocity_jet_error=max(self.velocity_jet_error,err)
        # Native inventory shifts are affine in the retained deep table.
        b=self.model.bulk;xi0=b.eos.xi.copy();ne=b.eos.gas(self.theta,self.eta)[6]
        scale=max(float(np.max(abs(b.eos.rinventory))),float(np.max(abs(b.eos.inventory[2])/ne)),1e-100);h=1e-5/scale
        values=[]
        for sign in [-1,1]:
            b.eos.xi=xi0+sign*h;cc=self.coefficients();values.append(cc)
        b.eos.xi=xi0
        dxi=(values[1]['emit']/self.scale-values[1]['loss']*X-values[0]['emit']/self.scale+values[0]['loss']*X)/(2*h)
        ratio=np.divide(values[1]['sc']-values[0]['sc'],2*h*c['sc'],out=np.zeros(self.n),where=c['sc']>0)
        ix=(dxi+ratio[:,None,None]*scat,dxi,ratio[None]*escaped)
        m=self.material;row=m.point(k);q=row['Q'];zero=np.zeros((4,self.n));field=np.zeros((5,self.n));maps=[]
        units=[np.maximum(q[0],1.),np.maximum(q[0]*C*C,1.),self.eu,self.nu]
        for j in range(4):
            z=zero.copy();z[j]=units[j];maps.append(primitive(m,k,z,field,c))
        field[0]=1/3;maps.append(primitive(m,k,zero,field,c))
        m.raw(k,zero,np.zeros_like(field),0.);bb=m.model.bulk;ut=bb.eos.gas(row['theta'],row['eta'])[2]
        imap={p:np.zeros(self.n) for p in ['dr','dv','dt','dy','dp']};imap['dt'][:self.nb]=bb.eos.inventory[1]/ut
        pg=bb.eos.gas(row['theta'],row['eta']);imap['dp'][:self.nb]=pg[4]*imap['dt'][:self.nb]-bb.eos.inventory[0]
        maps.append(imap)
        # Independent linearity of the direct variation on the actual saved state.
        test=primitive(m,k,self.motion[k],np.array([self.g['delta_u'][k],self.g['delta_log_lapse'][k],self.g['delta_lambda'][k],np.zeros(self.n),np.zeros(self.n)])/AMP,c)
        vol=(3*self.g['delta_u'][k]+self.g['delta_lambda'][k])/AMP
        for p in ['dr','dv','dt','dy']:
            estimate=sum(maps[j][p]*self.motion[k,j]/units[j] for j in range(4))+maps[4][p]*vol
            err=float(np.max(abs(test[p]-estimate))/max(np.max(abs(test[p])),1.));self.mapping_error=max(self.mapping_error,err);self.map_checks.append(dict(point=int(k),parameter=p,relative=err))
        def mapped(parts):return np.stack([sum(parts[p]*mp[p][:,None,None] for p in ['dr','dv','dt','dy']) for mp in maps],axis=-1)
        D=mapped(derivative);Db=mapped(dbound);De=np.stack([sum(desc[p]*mp[p][None] for p in ['dr','dv','dt','dy']) for mp in maps],axis=-1)
        D[...,5]+=ix[0];Db[...,5]+=ix[1];De[...,5]+=ix[2]
        row=dict(S=S,esc=esc,loss=c['loss'],sc=c['sc'],B=D[...,[2,3]],Bb=Db[...,[2,3]],Be=De[...,[2,3]],
            D=D[...,[0,1,4,5]],Db=Db[...,[0,1,4,5]],De=De[...,[0,1,4,5]],units=units[:2],em=em,coll=coll,bound=bound,escaped=escaped,
            pressure_map=np.stack([maps[j]['dp'] for j in [2,3]],axis=-1),pressure_drive=np.stack([maps[j]['dp'] for j in [0,1,4,5]],axis=-1),pressure=c['p'])
        self.cache[k]=row
        for j in list(self.cache):
            if j not in [k-1,k,k+1]:del self.cache[j]
        self.point_seconds+=time.monotonic()-start;self.point_count+=1
        return row

    def local(self,t):
        k=max(0,min(np.searchsorted(self.t,t,side='right')-1,15));f=(t-self.t[k])/(self.t[k+1]-self.t[k]);a=self.point(k);b=self.point(k+1)
        blend=lambda v:(1-f)*v[k]+f*v[k+1]
        volume=blend(3*self.g['delta_u']+self.g['delta_lambda'])/AMP;lapse=blend(self.g['delta_log_lapse'])/AMP
        motion=blend(self.motion);xi=np.r_[blend(self.xi),np.zeros(self.n-self.nb)]
        c={key:(1-f)*a[key]+f*b[key] for key in ['loss','sc','B','Bb','Be','pressure_map','S','esc']}
        c.update(q=np.zeros_like(self.I[0]),qb=np.zeros_like(self.I[0]),qe=np.zeros((3,self.n)),pressure_source=np.zeros(self.n))
        for w,r in [(1-f,a),(f,b)]:
            drive=np.stack([motion[0]/r['units'][0],motion[1]/r['units'][1],volume,xi],axis=-1)
            c['q']+=w*(np.einsum('nqfj,nj->nqf',r['D'],drive)+volume[:,None,None]*r['em']+lapse[:,None,None]*r['coll'])
            c['qb']+=w*(np.einsum('nqfj,nj->nqf',r['Db'],drive)+volume[:,None,None]*r['em']+lapse[:,None,None]*r['bound'])
            c['qe']+=w*(np.einsum('knj,nj->kn',r['De'],drive)+lapse[None]*r['escaped'])
            c['pressure_source']+=w*(np.einsum('nj,nj->n',r['pressure_drive'],drive)+volume*r['pressure'])
        c['mechanical']=np.diff(self.mechanical[k:k+2],axis=0)[0].T/(self.t[k+1]-self.t[k])/np.stack([self.eu,self.nu],axis=-1)
        return c

    def collision(self,c,x,g,source=False):
        p,q,e,b=super().collision(c,x,g,source)
        return p,q+c['mechanical'] if source else q,e,b

    def boundary_ports(self,t,x):
        k=max(0,min(np.searchsorted(self.t,t,side='left')-1,15));f=(t-self.t[k])/(self.t[k+1]-self.t[k]);I=(1-f)*self.I[k]+f*self.I[k+1]
        speed=((1-f)*self.g['delta_log_speed'][k]+f*self.g['delta_log_speed'][k+1])/AMP
        actual=x*self.scale+I*speed[:,None,None];ports=[]
        for j,area,mask in [(0,self.area[0],self.mu<0),(-1,self.area[-1],self.mu>0)]:
            values=4*np.pi*C*area*np.sum(actual[j,mask]*(self.w*self.mu)[mask,None]*self.num,axis=0)
            ports.append([float(values.sum()),float(values@self.E)])
        return np.array(ports)


# Retain the actual simultaneous SDIRK/Krylov owner. The added material input
# enters both stages and its own invariant ledger. Save collision-only transfer
# for the next real material sweep, plus both physical radial ports for GR.
runner=textwrap.dedent(inspect.getsource(mono.Response.run))
runner=replace(runner,'impulse=np.zeros(self.n)','impulse=np.zeros(self.n);transfer=np.zeros((self.n,2));ports=np.zeros((2,2));transfer_history=[];port_history=[]')
runner=replace(runner,"pressure=np.einsum('nj,nj->n',pmap,g)*self.volume","pressure=(np.einsum('nj,nj->n',pmap,g)+self.local(t)['pressure_source'])*self.volume\n        transfer_history.append(transfer.copy()*AMPLITUDE);port_history.append(ports.copy()*AMPLITUDE)")
runner=replace(runner,"g[:,0]*self.eu,g[:,1]*self.nu,impulse", "g[:,0]*self.eu+np.array([np.interp(t,self.t,v) for v in self.energy_offset.T]),g[:,1]*self.nu,impulse")
runner=replace(runner,"times=z['t'].tolist();records=list(z['moments']);", "transfer=z['collision_transfer'][-1]/AMPLITUDE;ports=z['radial_ports'][-1]/AMPLITUDE;transfer_history=list(z['collision_transfer']);port_history=list(z['radial_ports']);times=z['t'].tolist();records=list(z['moments']);")
runner=replace(runner,"qg=self.gas(c['q'],c['qb'],c['qe'])","qg=self.gas(c['q'],c['qb'],c['qe'])+c['mechanical']")
anchor="p2,g2,es2,_=self.collision(c,z,gz,True)"
runner=replace(runner,anchor,anchor+"\n        transfer+=h*((1-gamma)*(g1-c['mechanical'])+gamma*(g2-c['mechanical']))*np.stack([self.eu,self.nu],axis=-1)\n        ports+=h*((1-gamma)*self.boundary_ports(t+gamma*h,y)+gamma*self.boundary_ports(t+h,z))\n        ledger+=h*np.array([-np.sum(c['mechanical'][:,1]*self.nu),np.sum(c['mechanical'][:,0]*self.eu)])")
runner=replace(runner,"completed_steps=k+1,t=times,moments=records)","completed_steps=k+1,t=times,moments=records,transfer=transfer,ports=ports)")
runner=replace(runner,"ledger=ledger*AMPLITUDE,escape=escape*AMPLITUDE,radius_E=self.r)","ledger=ledger*AMPLITUDE,escape=escape*AMPLITUDE,radius_E=self.r,collision_transfer=transfer_history,radial_ports=port_history,energy_offset_reference=AMPLITUDE*self.energy_offset,energy_offset_t=self.t)")
runner=replace(runner,"row['passed']=max(error,species_error)<1e-8 and self.max_residual<1e-10","actual=g[:,0]*self.eu+np.array([np.interp(count*h,self.t,v) for v in self.energy_offset.T])\n    row.update(endpoint_material_nonrest_energy_erg=row['endpoint_material_reference_energy_erg'],endpoint_material_reference_energy_erg=float(np.sum(actual)*AMPLITUDE),endpoint_material_energy_L1_erg=float(np.sum(abs(actual))*AMPLITUDE),velocity_jet_relative=self.velocity_jet_error,primitive_mapping_relative=self.mapping_error,conserved_material_input_applied=True)\n    row['passed']=max(error,species_error)<1e-8 and self.max_residual<1e-10 and self.velocity_jet_error<.0001 and self.mapping_error<1e-10")
runner=replace(runner,'transfer_history=[];port_history=[]','transfer_history=[];port_history=[];photon_history=[];gas_history=[]')
runner=replace(runner,'transfer_history.append(transfer.copy()*AMPLITUDE);port_history.append(ports.copy()*AMPLITUDE)',
    'transfer_history.append(transfer.copy()*AMPLITUDE);port_history.append(ports.copy()*AMPLITUDE);photon_history.append(x*self.scale*AMPLITUDE);gas_history.append(g.copy()*AMPLITUDE)')
runner=replace(runner,"times=z['t'].tolist();records=list(z['moments']);", "photon_history=list(z['photon_history_scaled_occupation']) if 'photon_history_scaled_occupation' in z else [np.zeros_like(x),z['delta_packet_scaled_occupation']];gas_history=list(z['material_history']) if 'material_history' in z else [np.zeros_like(g),z['delta_material']];times=z['t'].tolist();assert len(times)==len(photon_history);records=list(z['moments']);")
runner=replace(runner,'transfer=transfer,ports=ports)','transfer=transfer,ports=ports,photon_history=photon_history,gas_history=gas_history)')
runner=replace(runner,'energy_offset_t=self.t)','energy_offset_t=self.t,photon_history_scaled_occupation=photon_history,material_history=gas_history)')
namespace=dict(vars(mono),OUT=OUT);exec(compile(runner,__file__,'exec'),namespace);Response.run=namespace['run']


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='d6ccddff9',
        claim='Return actual conserved baryon/momentum/inventory response to moving native collisions and simultaneously evolve photon packets with total material energy/H; save paired transfers and both radial energy ports for the next material/GR sweep.',
        decision='Measure how strongly the previously missing material response changes photons and the paired material force/energy. Then apply these transfers to material and GR, without calling one waveform sweep closed.',
        method='Reuse exact direct primitive variation. Local baryon/momentum and noncollisional material transport are interpolated from the accepted material history; energy/H are simultaneous unknowns with streaming and moving collisions in each SDIRK step. Account for volume once, without adding the old canonical-adiabatic thermal map twice.',
        velocity='Forward signed Richardson beta probes2e-7/1e-7 remain on the nonzero velocity remap branch. At exactly beta=0 the energy/H-to-velocity map vanishes; the initial prescribed velocity variation is also zero. Independently check derivative size sensitivity.',
        inventory='Keep native first-order advected inventory derivatives, including direct opacity/electron effects; do not treat composition work as pure heat.',
        reuse='Same531 cells,8 angles,152 frequencies,3.434ms,existing EOS/spectral banks and saved material/GR paths. No nonlinear background replay.',
        budgets=dict(check_seconds=40,pilot_seconds=65,production_seconds=900,CPU_threads=1,memory_GB=3,new_native_bank_calls=0),
        forecast='Prior monolithic three paths cost241.93s. This adds local conservative maps and signed velocity/inventory source jets at17 points. Measure point and per-step costs separately before dispatch;2x margin must fit900s. Late Krylov cost remains extrapolated.',
        gates=dict(owner=1e-10,mapping=1e-10,velocity_derivative=.0001,conservation=1e-8,linear=1e-10,time=.02,background=.02),
        stop='Only64/128 on reference128 and128 on reference64. Stop on gate or time cap; preserve source and failures. No automatic extra clock, horizon, physical grid, or full-closure claim.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(mono.__file__),Path(old.__file__),Path(matter.__file__),Path(primitive_owner.__file__),matter.OUT/'source-audit.json',mono.OUT/'result.json']}))
    (OUT/'expanded-run.py').write_text(runner);(OUT/'expanded-primitive.py').write_text(primitive_source)


def check():
    assert not (OUT/'check.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(40);m=Response(128);rows=[]
    rho,u,p,v,c2=sp.symbols('rho u p v c2',positive=True);W=1/sp.sqrt(1-v*v);H=rho*(c2+u)+p
    K=rho*c2*W*(W-1)+rho*u*W**2+p*(W**2-1)
    assert sp.simplify(K-(H*W**2-p-rho*c2*W))==0
    assert sp.simplify(sp.diff(K,v)-(rho*c2*v*W**3*(2*W-1)+2*(rho*u+p)*v*W**4))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Invertible Etilde=Eref-(a_ref-a_surface)*cx*c^2*B; exact special-relativistic nonrest energy K and velocity derivative. Fixed coefficients transfer the same baryon transport into the energy source; no heat or rest-mass contribution is dropped.'))
    for k in [0,8,16]:
        c=m.local(m.t[k]);x=np.zeros_like(m.I[0]);g=np.zeros((m.n,2));p,q,e,b=m.collision(c,x,g,True)
        paired=m.gas(p,b,e);energy=float(abs(np.sum(p*m.Eweight)+e[1].sum()+np.sum(paired[:,0]*m.eu))/max(np.sum(abs(p)*m.Eweight),np.sum(abs(paired[:,0])*m.eu),1.))
        species=float(abs(np.sum(b*m.Nweight)-np.sum(paired[:,1]*m.nu))/max(np.sum(abs(b)*m.Nweight),np.sum(abs(paired[:,1])*m.nu),1.))
        ports=m.boundary_ports(m.t[k],x);source=m.source(m.t[k])[1]/AMP
        # Radial-only balance removes spectral exits and metric frequency work.
        port_error=float(np.max(abs(ports[0]-ports[1]-np.array([source[3]+source[0],source[4]+source[1]-source[2]])))/max(np.max(abs(ports)),np.max(abs(source[:3])),1.))
        rows.append(dict(time=float(m.t[k]),energy=energy,species=species,radial_port=port_error))
    result=dict(classification='Counterexample candidate',passed=max(max(r['energy'],r['species'],r['radial_port']) for r in rows)<1e-10 and m.velocity_jet_error<.0001 and m.mapping_error<1e-10,
        rows=rows,velocity_jet_relative=m.velocity_jet_error,primitive_mapping_relative=m.mapping_error,mapping_checks=m.map_checks,point_seconds=m.point_seconds,seconds=time.monotonic()-start)
    write(OUT/'check.json',result);print(json.dumps(result));signal.alarm(0);assert result['passed']


def pilot():
    assert json.loads((OUT/'check.json').read_text())['passed'];assert not (OUT/'pilot.json').exists()
    start=time.monotonic();signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(65);rows=[]
    for steps in [64,128]:
        rows.append(Response(128).run(steps,f'pilot-{steps}',4))
        if not rows[-1]['passed']:break
    point=max(r['operator_point_seconds']/r['operator_points'] for r in rows)
    forecast=point*51+rows[0]['stepping_seconds']/4*60+rows[-1]['stepping_seconds']/4*252+25
    result=dict(classification='Counterexample candidate',rows=rows,forecast_seconds=forecast,upper_seconds=2*forecast,eligible=len(rows)==2 and all(r['passed'] for r in rows) and 2*forecast<900,seconds=time.monotonic()-start)
    write(OUT/'pilot.json',result);print(json.dumps({k:v for k,v in result.items() if k!='rows'}));signal.alarm(0)
    if result['eligible']:
        write(OUT/'execution-plan.json',dict(classification='Counterexample candidate',eligible=True,hard_cap_seconds=900,forecast_seconds=forecast,upper_seconds=2*forecast,
            paths=[[64,128,'pilot-64'],[128,128,'pilot-128'],[128,64,None]],
            correction='First map failed cancellation when a huge reference rest-mass transport offset was combined with tiny thermal energy. Use the exactly invertible nonrest energy coordinate and restore the offset in all physical readouts. Preserve the original check. Source pairing is checked before adding the independent mechanical source; port roundoff is normalized by the full frequency-work/source scale.',
            limits='One return sweep at prescribed material baryon/momentum/noncollision transport and original metric. No final coupled fixed point or GR result yet.',
            bindings={str(p):sha(p) for p in [Path(__file__),Path(mono.__file__),Path(old.__file__),Path(matter.__file__),Path(primitive_owner.__file__),OUT/'pilot.json',OUT/'check.json',matter.OUT/'source-audit.json']}))
        (OUT/'expanded-run-current.py').write_text(runner);(OUT/'expanded-primitive-current.py').write_text(primitive_source)


def production():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'execution-plan.json').read_text());assert plan['eligible']
    for p,h in plan['bindings'].items():assert sha(p)==h,p
    start=time.monotonic();signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(plan['hard_cap_seconds']);rows=[]
    try:
        for steps,ref,restart in plan['paths']:
            rows.append(Response(ref).run(steps,f'steps-{steps}-reference-{ref}',restart=restart))
            if not rows[-1]['passed']:break
        comparisons={};indices=[0,1,2,3,5,6]
        if len(rows)==3 and all(r['passed'] for r in rows):
            def history(steps,ref):
                d=np.load(OUT/f'steps-{steps}-reference-{ref}.npz');t=np.linspace(0,flow.old.END,17);ids=[int(np.argmin(abs(d['t']-s))) for s in t];assert np.max(abs(d['t'][ids]-t))<1e-18
                return d['moments'][ids][:,indices]
            fine=history(128,128);norm=np.maximum(np.max(np.sum(abs(fine),axis=2),axis=0),1.)
            for name,v in [('time',history(64,128)),('background',history(128,64))]:comparisons[name]=(np.max(np.sum(abs(v-fine),axis=2),axis=0)/norm).tolist()
        result=dict(classification='Counterexample candidate',passed=len(rows)==3 and all(r['passed'] for r in rows) and max((max(v) for v in comparisons.values()),default=1)<.02,
            comparisons=comparisons,comparison_order=['photon_energy','material_reference_energy','neutral_H','collision_momentum_impulse','photon_radial_pressure','material_pressure'],paths=rows,seconds=time.monotonic()-start,
            actual_material_response_returned_to_photons=True,returned_collision_transfer_applied_to_material=False,returned_source_applied_to_GR=False,coupled_fixed_point_verified=False,final_charge_solved=False)
        write(OUT/'result.json',result);print(json.dumps(result),flush=True)
    except Exception as exc:
        write(OUT/'production-failure.json',dict(error=repr(exc),completed_paths=rows,seconds=time.monotonic()-start));raise
    finally:signal.alarm(0)


def history_receipt():
    assert not (OUT/'result.json').exists();p=json.loads((OUT/'execution-plan.json').read_text());assert p['eligible']
    p['bindings'][str(Path(__file__))]=sha(__file__)
    p['history_receipt']='Save full photon and material histories at existing canonical17 times, plus the pilot endpoint. This is required to evaluate the next coupled-equation residual without replaying an accepted path. The pilot has exactly the initial and current endpoint; reconstruct those from its saved exact arrays. Stage equations, sources, Krylov gates and grids are unchanged. Added RAM below120MB per path; compressed artifacts retained for reuse.'
    write(OUT/'execution-plan.json',p);(OUT/'expanded-run-current.py').write_text(runner)


if __name__=='__main__':globals()[sys.argv[1]]()
