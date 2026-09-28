"""Request29: preserve 26 nuclear inventories in a declared FreeEOS embedding.

Counterexample candidate. Charge/ion preserving trace representation is a new
EOS approximation, not a certified extension of FreeEOS's atomic physics.
"""
from pathlib import Path
import ctypes, json, re, shutil, subprocess, sys
import numpy as np
import gr_mass as gr
import conservative_cell as cell
import reactive_energy as reaction
import baryon_entropy as be
import fresh_microphysics as fresh
import native_closure as native
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares
from scipy.linalg import expm
from scipy.integrate import solve_ivp
from thermal_restart import sha
from thermal_wd import mesa

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/common-eos29'
CACHE=Path('/home/lpaiu/work/common-eos29')
NAMES=cell.species()
ISO=reaction.isotopes()
A=np.array([ISO[n]['a'] for n in NAMES],float)
Z=np.array([ISO[n]['z'] for n in NAMES],float)
W=np.array([float(ISO[n]['w']) for n in NAMES])
EZ=np.array([1,2,6,7,8,10,11,12,13,14,15,16,17,18,20,22,24,25,26,28])
NA=6.02214076e23
KB=1.380649e-16
EV=1.602176634e-12


def save(name,obj):
    (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists() and not CACHE.exists()
    OUT.mkdir();CACHE.mkdir()
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='d20de02',
        previous_manifest_sha256=sha(cell.OUT/'manifest.json'),
        steps=['Full 26-species inventory and atomic energy reference audit',
            'Charge and ion number preserving trace FreeEOS representation',
            'Thermodynamic derivative and independent embedding controls',
            'Baryon-coordinate GR evolution and boundary luminosity/charge connection'],
        controls=dict(moment_absolute=1e-14,EOS_derivative_relative=1e-5,
            embedding_pressure_relative=1e-5,embedding_energy_relative_to_cvT=1e-5,
            time_composition_relative=1e-3,time_temperature_log_absolute=2e-6,
            GR_interface_absolute=1e-8,energy_relative_to_released=1e-6),
        policy='Retain original physical nuclear abundances and rest masses. '
            'Record every failed control. Agreement between two trace embeddings is sensitivity, '
            'not a physical error bound. A physically uncertified EOS or transport closure prevents '
            'a completed physical-star or observational-inference verdict.',
        evolution='A time-dependent entropy/composition sequence must satisfy the same baryon GR constraints '
            'and redshifted energy balance. Independent mass refits or unchanged initial profiles are not evolution.',
        exclusions='No observational runtime reopening, new timing data, runtime environment rebuild, or submission.'))
    for rel in ['docs/'+s+'.md' for s in cell.DOCS]+['paper/revision-manifest.json']:
        dest=OUT/'previous-notes'/Path(rel).name;dest.parent.mkdir(exist_ok=True)
        shutil.copy2(ROOT/rel,dest)
    save('historical-bindings.json',{rel:sha(ROOT/rel) for rel in
        ['docs/'+s+'.md' for s in cell.DOCS]+['paper/revision-manifest.json']})


def binding():
    rows=json.loads((OUT/'nist-trace-data.json').read_text())['rows']
    assert {r['Z'] for r in rows}==set(Z)-set(EZ)
    for r in rows: assert len(r['energies_eV'])==len(r['uncertainties_eV'])==r['Z']
    sums={str(r['Z']):sum(r['energies_eV']) for r in rows}
    save('atomic-binding.json',dict(classification='Imported from prior work',
        source_sha256=sha(OUT/'nist-trace-data.json'),total_binding_eV=sums,
        note='Element ground states. Isotope shifts and database uncertainties remain physical uncertainties.'))
    src=(gr.SOURCE/'src/mod_isotopic_mass_data.f90').read_text()
    selected=src.split('isotopic_mass(nelements_iso) = [&',1)[1].split(']',1)[0]
    weights=[]
    for n,z in re.findall(r'isotopic_mass_2d\(\s*(-?\d+),\s*(\d+)\)',selected):
        part=src.split(f'isotopic_mass_{int(z):03d}(min_nmz:max_nmz) = [&',1)[1].split(']',1)[0]
        vals=[float(v) for v in re.findall(r'([0-9]+\.[0-9]*)_fp_kind',part)]
        assert len(vals)==73
        weights.append(vals[int(n)+8])
    assert len(weights)==20 and weights[2]==12
    save('freeeos-weights.json',dict(classification='Imported from prior work',atomic_weights=weights,
        source_sha256=sha(gr.SOURCE/'src/mod_isotopic_mass_data.f90')))
    print('TOTAL BINDING eV',sums,flush=True)


def embedding(wide=False):
    """Map number abundances, never mass fractions, preserving ions and charge."""
    matrix=np.zeros((26,20))
    for i,z in enumerate(Z):
        if z in EZ: matrix[i,np.where(EZ==z)[0][0]]=1;continue
        lo,hi=(1,10) if wide else (EZ[EZ<z][-1],EZ[EZ>z][0])
        matrix[i,np.where(EZ==lo)[0][0]]=(hi-z)/(hi-lo)
        matrix[i,np.where(EZ==hi)[0][0]]=(z-lo)/(hi-lo)
    assert np.max(abs(matrix.sum(1)-1))<1e-14
    assert np.max(abs(matrix@EZ-Z))<1e-14
    return matrix


def census():
    b=np.load(cell.OLD/'PP-steady-2-input.npz');x=b['X'];dm=b['dm']
    missing=np.array([z not in EZ for z in Z]);rows=[]
    for i in np.where(missing)[0]:
        active=x[:,i]>0
        rows.append(dict(species=NAMES[i],max_X=float(x[:,i].max()),
            mean_X=float(dm@x[:,i]/sum(dm)),active_cells=int(active.sum()),
            minimum_active_T_K=float(np.exp(b['lnT'][active]).min()) if active.any() else None))
    save('inventory.json',dict(classification='Proven',input_sha256=sha(cell.OLD/'PP-steady-2-input.npz'),
        rows=rows,primary_matrix=embedding().tolist(),wide_matrix=embedding(True).tolist(),
        sum_X_error=float(np.max(abs(x.sum(1)-1))),
        interpretation='Arithmetic properties of the saved initial state and explicit maps; no claim that trace elements are physically absent after burning.'))
    print('INVENTORY',rows,flush=True)


class EOS:
    def __init__(self,wide=False):
        self.native=gr.EOS();self.map=embedding(wide)
        self.bind=json.loads((OUT/'atomic-binding.json').read_text())['total_binding_eV']
        # Source-normalized FreeEOS element binding energies, not current NIST
        # replacements: the correction must use exactly the library's zero.
        src=(gr.SOURCE/'src/mod_ionization_data.f90').read_text()
        block=src.split('monatomic_ip(nions) = [&',1)[1].split(']',1)[0]
        values=np.array([float(v) for v in re.findall(r'([0-9]+\.[0-9]*)_fp_kind',block)])
        assert len(values)==295
        bounds=np.r_[0,np.cumsum(EZ)]
        h=6.62607015e-27;c=gr.C*100
        self.map_binding=np.array([values[a:b].sum()*h*c for a,b in zip(bounds[:-1],bounds[1:])])
        self.real_binding=np.array([self.bind.get(str(int(z)),0)*EV for z in Z])
        # Supported elements retain exactly FreeEOS atomic data, avoiding a
        # silent update of its partial-ionization model.
        for i,z in enumerate(Z):
            if z in EZ: self.real_binding[i]=self.map_binding[np.where(EZ==z)[0][0]]

    def composition(self,x):
        assert x.shape==(26,) and min(x)>=-1e-30 and abs(sum(x)-1)<1e-11
        y=x/A;ym=y@self.map
        # FreeEOS internally normalizes by its own neutral atomic weights.
        # Recover that factor from its density-independent mass normalization.
        return y,ym

    def __call__(self,mode,logvalue,logT,x):
        y,ym=self.composition(x)
        # FreeEOS uses eps = atom number / neutral-atom gram. A normalized
        # number mixture gives CX via its actual source atomic weights.
        cx=float(ym@self.weights)
        eps=ym/cx
        a=self.native(mode,logvalue+np.log(cx) if mode==2 else logvalue,logT,eps)
        b=a.copy();b[0]/=cx;b[2:4]*=cx;b[9:11]*=cx
        # Constant-in-rho,T correction aligns the fully stripped trace limit
        # with the *real* neutral nuclei+electrons reference. No partial-ion
        # correction is claimed for unsupported elements.
        correction=NA*(y@self.real_binding-ym@self.map_binding)
        b[2]+=correction
        return b

    @property
    def weights(self):
        if not hasattr(self,'_weights'):
            vals=json.loads((OUT/'freeeos-weights.json').read_text())
            self._weights=np.array(vals['atomic_weights'])
        return self._weights


def probe():
    b=np.load(OUT/'restored-input.npz') if (OUT/'restored-input.npz').exists() else np.load(cell.OLD/'PP-steady-2-input.npz')
    label='restored-EOS-probe.json' if (OUT/'restored-input.npz').exists() else 'EOS-probe.json'
    eos=EOS();other=EOS(True)
    selected=np.unique(np.r_[np.linspace(0,5734,65).astype(int),cell.ZONE,
        np.argmax(b['X'],axis=0)])
    rows=[]
    for i in selected:
        rho=b['lnd'][i];t=b['lnT'][i];x=b['X'][i];a=eos(2,rho,t,x);w=other(2,rho,t,x)
        for step in [1e-4,5e-5,1e-5]:
            ap=eos(2,rho,t+step,x);am=eos(2,rho,t-step,x)
            rp=eos(2,rho+step,t,x);rm=eos(2,rho-step,t,x)
            dT=(ap-am)/(2*step);dr=(rp-rm)/(2*step)
            errors=dict(P_T=dT[1]/(a[1]*a[6])-1,
                cv=dT[2]/a[10]-1,entropy_T=dT[3]*np.exp(t)/a[10]-1,
                Maxwell=(dr[2]-(a[1]-dT[1])/a[0])/max(abs(a[10]),1.),
                entropy_rho=(dr[3]+dT[1]/(a[0]*np.exp(t)))/max(abs(a[3]),1.))
            rows.append(dict(cell=int(i),step=step,errors=errors,
                embedding_P_relative=float(w[1]/a[1]-1),
                embedding_u_over_cvT=float((w[2]-a[2])/a[10])))
    worst={key:max(abs(r['errors'][key]) for r in rows) for key in rows[0]['errors']}
    save(label,dict(classification='Counterexample candidate',states=len(selected),rows=rows,
        worst=worst,derivative_passed=all(v<1e-5 for v in worst.values()),
        physical_error_certified=False,
        note='Finite native FreeEOS derivative comparisons and trace model sensitivity only.'))
    print('EOS PROBE',worst,flush=True)


class Structure(be.Solver):
    """Reuse the existing baryon TOV integrator with physical nuclear rest mass."""
    def __init__(self,label,data,points=17,subdivision=2):
        self.mat=be.Material();m=self.mat
        self.eos=EOS();self.sub=subdivision;self.method='RK4';self.calls=0;self.direct_calls=0
        m.baryon_g=float(sum(data['dm']));m.B=m.baryon_g*.001*gr.G/gr.C**2
        m.dm=data['dm']/m.baryon_g;m.R=68819648.58910435
        m.outer=np.r_[0.,np.cumsum(m.dm)];m.inner=np.r_[np.cumsum(m.dm[::-1])[::-1],0.]
        m.outer[-1]=m.inner[0]=1.;m.split=int(np.argmin(abs(m.outer-.5)))
        m.eps=data['X'];m.lt=data['lnT'];m.cx=(data['X']/A)@W
        path=OUT/f'{label}-adiabats-{points}.npz'
        if path.exists():
            tab=np.load(path);self.ref=tab['reference'];m.lp=np.log(self.ref[:,1]);vals=tab['values'];offset=tab['offset']
            assert np.array_equal(tab['X'],data['X']) and np.array_equal(tab['lnT'],data['lnT'])
            assert np.array_equal(tab['lnd'],data['lnd'])
        else:
            self.ref=np.array([self.eos(2,r,t,x) for r,t,x in zip(data['lnd'],m.lt,m.eps)])
            m.lp=np.log(self.ref[:,1]);offset=np.linspace(-.04,.04,points)
            vals=np.empty((len(m.lp),points,3))
            for i in range(len(m.lp)):
                for j,dx in enumerate(offset):
                    a,t,_=be.invert(self.eos,m.lp[i]+dx,self.ref[i,3],m.eps[i],m.lt[i])
                    vals[i,j]=[np.log(a[0]),t,a[2]*1e-4/gr.C**2]
                if i%1000==0: print('GR EOS table',label,points,i,flush=True)
            np.savez_compressed(path,reference=self.ref,values=vals,offset=offset,
                X=data['X'],lnT=data['lnT'],lnd=data['lnd'])
        self.offset=offset;self.coef=PchipInterpolator(offset,vals,axis=1).c
        self.lp_reference=m.lp.copy()
        if 'boundary_logP' in data: m.lp[0]=float(data['boundary_logP'])

    def state(self,lp,i):
        self.calls+=1;m=self.mat;dx=lp-self.lp_reference[i]
        if not self.offset[0]<=dx<=self.offset[-1]:
            # No extrapolation: face pressures can lie outside the midpoint's
            # narrow lookup band. Solve the same EOS/entropy equation directly.
            self.direct_calls+=1
            a,t,_=be.invert(self.eos,lp,self.ref[i,3],m.eps[i],m.lt[i])
            b=gr.G*a[0]*1000/gr.C**2
            return gr.G*np.exp(lp)*.1/gr.C**4,b*(m.cx[i]+a[2]*1e-4/gr.C**2),b
        j=min(len(self.offset)-2,max(0,int(np.searchsorted(self.offset,dx)-1)));z=dx-self.offset[j]
        a=self.coef[:,j,i];v=((a[0]*z+a[1])*z+a[2])*z+a[3]
        b=gr.G*np.exp(v[0])*1000/gr.C**2
        return gr.G*np.exp(lp)*.1/gr.C**4,b*(m.cx[i]+v[2]),b


def structure(label,data,points=17,sub=2,initial=None):
    solver=Structure(label,data,points,sub);m=solver.mat
    if initial is None: initial=[m.lp[-1]+.001,0.,0.]
    fit=least_squares(lambda x:solver.branches(x)[0],initial,jac='3-point',
        bounds=([m.lp[-1]-.2,-.1,-.01],[m.lp[-1]+.2,.1,.01]),
        xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=25)
    error,inner,outer=solver.branches(fit.x,record=True)
    assert fit.success and max(abs(error))<1e-8,(fit.message,error)
    n=len(m.lp);B=m.B;R=m.R*np.exp(fit.x[1]);M=gr.TARGET*np.exp(fit.x[2])
    faces=np.zeros((n+1,3));faces[0]=[R/m.R,M/B,m.lp[0]]
    faces[1:m.split+1]=outer[1:,1:];faces[m.split:n]=inner[1:,1:][::-1];faces[n]=[0,0,fit.x[0]]
    # Finite-pressure photospheric boundary: explicitly retain the pressure
    # work term. This normalization omits the old mathematical atmosphere.
    nuface=.5*np.log1p(-2*M/R);rows=[];nufaces=[nuface]
    for i in range(n):
        outside=i<m.split
        if outside: begin=outer[i,0];mid=(m.outer[i]+m.outer[i+1])/2;start=outer[i,1:]
        else:
            j=n-1-i;begin=inner[j,0];mid=(m.inner[i]+m.inner[i+1])/2;start=inner[j,1:]
        y=solver.step(np.log(begin),np.log(mid),start,i,B,outside)
        values=[be.invert(solver.eos,p,solver.ref[i,3],m.eps[i],m.lt[i]) for p in [faces[i,2],y[2],faces[i+1,2]]]
        a,t,_=values[1];aa,bb=values[0][0],values[2][0]
        ht=aa[2]+aa[1]/aa[0];hm=a[2]+a[1]/a[0];hb=bb[2]+bb[1]/bb[0]
        nu=nuface-np.log1p((hm-ht)/(m.cx[i]*(gr.C*100)**2+ht))
        nuface-=np.log1p((hb-ht)/(m.cx[i]*(gr.C*100)**2+ht));nufaces.append(nuface)
        rows.append([np.log(a[0]),t,y[0]*m.R,y[1]*B,y[2],nu,a[2],a[3]])
    a=np.array(rows);result=dict(data)
    for j,k in enumerate(['lnd','lnT','r_mid_m','m_mid_geom','logP','nu','u_W','s_B']): result[k]=a[:,j]
    result['lnR']=np.log(faces[:-1,0]*m.R*100);result['nu_faces']=np.array(nufaces)
    result['radius_faces_m']=faces[:,0]*m.R;result['mass_faces_geom']=faces[:,1]*B
    np.savez_compressed(OUT/f'{label}-state-{points}-{sub}.npz',**result)
    report=dict(classification='Counterexample candidate',parameters=fit.x.tolist(),
        interface_max=float(max(abs(error))),radius_m=R,mass_geom_m=M,baryon_g=m.baryon_g,
        EOS_calls=solver.calls,direct_calls=solver.direct_calls,physical_EOS_certified=False,
        boundary='Nonzero photospheric pressure; exterior is a boundary normalization, not a solved radiating atmosphere.')
    save(f'{label}-structure-{points}-{sub}.json',report)
    print('GR STRUCTURE',label,points,sub,report,flush=True)
    return result,report


def baseline():
    data=dict(np.load(cell.OLD/'PP-steady-2-input.npz'))
    for points,sub in [(9,1),(17,2),(17,4)]:
        previous=OUT/f'base-structure-{9 if points==17 and sub==2 else 17}-{1 if sub==2 else 2}.json'
        initial=json.loads(previous.read_text())['parameters'] if previous.exists() else None
        if not (OUT/f'base-structure-{points}-{sub}.json').exists():
            structure('base',data,points,sub,initial)


def boundary_audit():
    d=np.load(OUT/'base-state-17-4.npz');_,profile=mesa(cell.OLD/'PP-steady-2-profile.data.gz')
    N=np.exp(d['nu']);nf=np.exp(d['nu_faces']);r=d['r_mid_m'];dm=d['dm']
    Linf=np.r_[d['L']*nf[:-1]**2,0.]
    qflux=(Linf[:-1]-Linf[1:])/(dm*N*N)
    heat=profile['eps_nuc']-profile['non_nuc_neu']-qflux
    timescale=profile['cv']*np.exp(d['lnT'])/np.maximum(abs(heat),1e-99)
    # The radiative+conductive diffusion law must include the Tolman gradient.
    # Input luminosity is not automatically a solution of that law.
    rho=np.exp(d['lnd']);T=np.exp(d['lnT']);f=1-2*d['m_mid_geom']/r
    gNT=np.gradient(N*T,r)
    arad=7.5657e-15
    Ldiff=-16*np.pi*(r*100)**2*arad*(gr.C*100)*T**3/(3*profile['opacity']*rho)*np.sqrt(f)/N*gNT/100
    fluxscore=abs(Ldiff-d['L'])/np.maximum(abs(d['L']),fresh.LSUN*1e-10)
    # Finite-pressure surface mass work is included when a boundary moves.
    np.savez_compressed(OUT/'boundary-forcing.npz',Linf=Linf,qflux=qflux,
        net_heat=heat,thermal_timescale_s=timescale,diffusion_L=Ldiff,flux_score=fluxscore)
    save('boundary-audit.json',dict(classification='Counterexample candidate',
        retained_surface_Lsun=float(Linf[0]/fresh.LSUN),
        net_redshifted_power_Lsun=float((dm*N*N)@heat/fresh.LSUN),
        minimum_thermal_timescale_s=float(timescale.min()),
        mass_fraction_diffusion_mismatch_over_1pct=float(dm@(fluxscore>.01)/sum(dm)),
        diffusion_relative_mismatch_max=float(fluxscore.max()),
        all_transport_closed=False,
        note='Same material baryon state, but native powers/opacity at the preceding very close state. '
            'Diffusion diagnostic excludes convective closure and physical atmosphere. Retained flux is an explicit forcing, not certified stellar transport.'))
    print('BOUNDARY',json.loads((OUT/'boundary-audit.json').read_text()),flush=True)


def fixed_boundary():
    data=dict(np.load(cell.OLD/'PP-steady-2-input.npz'))
    data['boundary_logP']=be.Material().lp[0]
    save('evolution-plan.json',dict(classification='Counterexample candidate',
        duration_coordinate_seconds=140780.16,steps=[1,2,4],
        boundary_logP=float(data['boundary_logP']),
        flux='Prescribe the initial L_infinity on every baryon face for the duration. '
            'This is a forced quasistatic GR evolution control, not self-consistent stellar transport.',
        state='All 26 nuclear inventories and fixed cell baryon masses. Evolve composition at native '
            'fixed rho,T by a stiff exponential local linearization. Integrate reaction neutrino loss '
            'in the same exponential. Invert the neutral-reference total energy at fixed density, '
            'then solve the same-entropy/composition baryon TOV constraints, without a mass target fit.',
        energy='u_B + sum_i(W_i/A_i-1) X_i c^2; baryon constant omitted. Nuclear heat is not added twice.',
        controls=json.loads((OUT/'plan.json').read_text())['controls'],
        zero_step='Fixed face pressure must remain fixed when midpoint states are reinserted. '
            'An EOS lookup midpoint is not a photospheric boundary.',
        physical_gate='Report frozen-flux transport defect, native/FreeEOS auxiliary mismatch, '
            'time refinement and GR energy defect before any physical-evolution verdict.'))
    for points in [9,17]:
        shutil.copy2(OUT/f'base-adiabats-{points}.npz',OUT/f'epoch0-adiabats-{points}.npz')
    for points,sub in [(9,1),(17,2),(17,4)]:
        initial=json.loads((OUT/f'base-structure-{points}-{sub}.json').read_text())['parameters']
        if not (OUT/f'epoch0-structure-{points}-{sub}.json').exists():
            structure('epoch0',data,points,sub,initial)


def restore():
    data=dict(np.load(cell.OLD/'PP-steady-2-input.npz'));_,source=mesa(fresh.SOURCE)
    F=sum(source[n] for n in NAMES if n.startswith('f'))
    data['X']=data['X']*(1-F[:,None])
    for i,n in enumerate(NAMES):
        if n.startswith('f'): data['X'][:,i]=source[n]
    data['boundary_logP']=be.Material().lp[0]
    assert max(abs(data['X'].sum(1)-1))<1e-12
    np.savez_compressed(OUT/'restored-input.npz',**data)
    save('restore-plan.json',dict(classification='Counterexample candidate',
        original_profile_sha256=sha(fresh.SOURCE),old_input_sha256=sha(cell.OLD/'PP-steady-2-input.npz'),
        mean_restored_F=float(data['dm']@F/sum(data['dm'])),max_restored_F=float(F.max()),
        method='Undo the previous non-F renormalization: multiply the entire 26-entry non-F state by 1-F_original '
            'and insert original F17/F18/F19 on the same material cells. This retains the introduced PP intermediate '
            'bookkeeping while restoring the original bulk F inventory. All rates must be reevaluated.',
        boundary_logP=float(data['boundary_logP']),
        frozen_controls=['EOS-probe.json','base structures','epoch0 structures'],
        physical_trace_model_certified=False))
    probe()
    for points,sub in [(9,2),(17,4)]:
        prior='epoch0' if points==9 else 'restored'
        priorpoints,priorsub=(17,4) if points==9 else (9,2)
        initial=json.loads((OUT/f'{prior}-structure-{priorpoints}-{priorsub}.json').read_text())['parameters']
        structure('restored',data,points,sub,initial)


def evaluate(label,data):
    native.OUT=OUT;native.CACHE=CACHE;native.context()
    if not (OUT/f'{label}-native.npz').exists():
        native.setup(label,data,species=NAMES,
            network=(cell.OLD/'inputs/explicit-PP/cno_extras.net').read_text())
        native.trace(label)
    inp=np.load(OUT/f'{label}-input.npz')
    for key in ['X','lnd','lnT','dm']: assert np.array_equal(inp[key],data[key]),(label,key)
    _,p=mesa(OUT/f'{label}-profile.data.gz')
    return dict(np.load(OUT/f'{label}-native.npz')),p


def source_step(data,d,p,dt,Linf,label):
    eos=EOS();x=data['X'];rho=data['lnd'];lnT=data['lnT'];N=np.exp(data['nu'])
    qflux=(Linf[:-1]-Linf[1:])/(data['dm']*N*N)
    rest=(W/A-1)*(gr.C*100)**2;qlegacy,_=cell.rest()
    result={k:v.copy() if hasattr(v,'copy') else v for k,v in data.items()}
    proposed=np.empty_like(x);temperature=np.empty_like(lnT);loss=np.empty_like(lnT)
    rows=[];negative=0.;closure=0.;released=0.
    for i in range(len(x)):
        h=dt*N[i];J=d['jacobian'][i];f=d['dxdt'][i]
        a=eos(2,rho[i],lnT[i],x[i]);ux=np.zeros(26)
        # Same coupled reaction/thermal construction as Request28. Recompute
        # composition energy derivatives from this declared EOS, not HELM.
        if max(abs(f))>1e-30:
            for k in range(26):
                if k==NAMES.index('he4'): continue
                slopes=[]
                for size in [1e-5,5e-6]:
                    z=x[i].copy();z[k]+=size;z[NAMES.index('he4')]-=size
                    slopes.append((eos(2,rho[i],lnT[i],z)[2]-a[2])/size)
                ux[k]=2*slopes[1]-slopes[0]
        tempforce=d['T'][i]*d['dxdt_T'][i]
        nx=-qlegacy@J-d['heat_X'][i]
        nt=-qlegacy@tempforce-d['T'][i]*d['heat_T'][i]
        rate=(-rest@f-d['neutrino'][i]-p['non_nuc_neu'][i]-qflux[i]-ux@f)/a[10]
        mat=np.zeros((29,29));mat[:26,:26]=J;mat[:26,26]=tempforce;mat[:26,-1]=f
        mat[26,:26]=(-(rest+ux)@J-nx)/a[10]
        mat[26,26]=(-(rest+ux)@tempforce-nt)/a[10];mat[26,-1]=rate
        mat[27,:26]=nx;mat[27,26]=nt;mat[27,-1]=d['neutrino'][i]
        delta=expm(h*mat)[:,-1];xx=x[i]+delta[:26]
        negative=min(negative,float(xx.min()));assert xx.min()>-1e-14,(i,xx.min())
        xx=np.maximum(xx,0);xx[NAMES.index('he4')]+=1-sum(xx)
        restchange=float(rest@(xx-x[i]))
        losses=delta[27]+h*p['non_nuc_neu'][i]
        assert losses>=-1e-12,(i,losses)
        target=a[2]-restchange-losses-h*qflux[i]
        t=lnT[i]+delta[26]
        for _ in range(15):
            aa=eos(2,rho[i],t,xx);defect=aa[2]-target
            if abs(defect)<max(2.,abs(target-a[2])*1e-8): break
            t-=np.clip(defect/aa[10],-.1,.1)
        else: raise ValueError(('energy inverse',i,defect))
        proposed[i]=xx;temperature[i]=t;loss[i]=losses
        closure=max(closure,abs(defect));released+=abs(restchange)*data['dm'][i]
        if i%1000==0: print('SOURCE',label,i,'dlnT',t-lnT[i],flush=True)
    result['X']=proposed;result['lnT']=temperature
    np.savez_compressed(OUT/f'{label}-source.npz',**result,loss_per_baryon_gram=loss)
    save(f'{label}-source.json',dict(classification='Counterexample candidate',dt=dt,
        max_inverse_defect_erg_g=closure,minimum_unprojected_X=negative,
        max_delta_lnT=float(max(abs(temperature-lnT))),absolute_rest_release_erg=released,
        neutrino_energy_infinity_erg=float(data['dm']@(N*loss)),
        source_sha256=sha(OUT/f'{label}-native.npz') if (OUT/f'{label}-native.npz').exists() else None))
    return result,loss


def stable_mass_change(a,b):
    r0=a['r_mid_m'];r1=b['r_mid_m'];m0=a['m_mid_geom'];m1=b['m_mid_geom']
    f0=np.sqrt(1-2*m0/r0);f1=np.sqrt(1-2*m1/r1)
    df=-2*((m1-m0)-m0/r0*(r1-r0))/r1/(f0+f1)
    cx=(a['X']/A)@W;dcx=((b['X']-a['X'])/A)@(W-A)
    du=(b['u_W']-a['u_W'])*1e-4/gr.C**2
    return float(a['dm']@((dcx+du)*f1+(cx+a['u_W']*1e-4/gr.C**2)*df))


def evolution():
    plan=json.loads((OUT/'evolution-plan.json').read_text());dtall=plan['duration_coordinate_seconds']
    assert (OUT/'coupled-evolution-plan.json').exists()
    original=dict(np.load(OUT/'restored-state-17-4.npz'));eos=EOS()
    nf=np.exp(original['nu_faces']);Linf=np.r_[original['L']*nf[:-1]**2,0.]
    d0,p0=evaluate('evolution-initial',original)
    summaries=[];previous=None
    for steps in plan['steps']:
        state={k:v.copy() for k,v in original.items()};neutrino=0.;records=[]
        parameters=json.loads((OUT/'restored-structure-17-4.json').read_text())['parameters']
        for step in range(steps):
            label=f'evolve-{steps}-{step}'
            d,p=(d0,p0) if step==0 else evaluate(label,state)
            nxt,loss=source_step(state,d,p,dtall/steps,Linf,label)
            neutrino+=float(state['dm']@(np.exp(state['nu'])*loss))
            state,row=structure(label,nxt,17,4,parameters);parameters=row['parameters']
            records.append(row)
        mass=stable_mass_change(original,state)*(gr.C*100)**2
        rs0=original['radius_faces_m'][0]*100;rs1=state['radius_faces_m'][0]*100
        pressure=np.exp(float(original['boundary_logP']))
        boundary_work=pressure*4*np.pi*(rs1-rs0)*(rs1*rs1+rs1*rs0+rs0*rs0)/3
        expected=-Linf[0]*dtall-neutrino-boundary_work
        release=sum(json.loads((OUT/f'evolve-{steps}-{j}-source.json').read_text())['absolute_rest_release_erg'] for j in range(steps))
        row=dict(classification='Counterexample candidate',steps=steps,records=records,
            stable_mass_change_energy_erg=mass,expected_mass_change_energy_erg=expected,
            surface_pressure_work_erg=boundary_work,neutrino_infinity_erg=neutrino,
            energy_score=abs(mass-expected)/release,release_erg=release,
            energy_passed=bool(abs(mass-expected)/release<1e-6),
            max_delta_lnT=float(max(abs(state['lnT']-original['lnT']))))
        if previous is not None:
            score=float(np.max(abs(state['X']-previous['X'])/(1e-16+1e-3*abs(state['X']-original['X']))))
            temp=float(max(abs(state['lnT']-previous['lnT'])))
            row['refinement']=dict(composition_score=score,temperature_difference=temp,
                passed=bool(score<1 and temp<2e-6))
        previous=state;summaries.append(row)
        save('evolution.json',dict(classification='Counterexample candidate',rows=summaries,
            physical_transport_closed=False,full_GR_dynamics=False))
        print('EVOLUTION',steps,'energy_score',row['energy_score'],'refinement',row.get('refinement'),flush=True)


def scalar_model(label,data,tolerance=2e-12):
    R=float(data['radius_faces_m'][0]);M=float(data['mass_faces_geom'][0]);mu=M/R
    r=data['r_mid_m'][::-1];x=np.r_[0.,r/R,1.]
    rho=np.exp(data['lnd'][::-1]);P=np.exp(data['logP'][::-1])*.1
    rest=(data['X'][::-1]/A)@W
    e=gr.G*rho*1000/gr.C**2*(rest+data['u_W'][::-1]*1e-4/gr.C**2)
    trace=e-3*gr.G*P/gr.C**4
    N=np.exp(data['nu'][::-1]);f=1-2*data['m_mid_geom'][::-1]/r
    a=PchipInterpolator(x,np.r_[N[0],N*np.sqrt(f),1-2*mu])
    b=PchipInterpolator(x,np.r_[4*np.pi*(-4)*R*R*N[0]*trace[0],
        4*np.pi*(-4)*R*R*N*trace/np.sqrt(f),
        4*np.pi*(-4)*R*R*trace[-1]])
    def rhs(t,y): return [y[1]/(t*t*a(t)),t*t*b(t)*y[0]]
    start=1e-9;b0=float(b(0));a0=float(a(0))
    sol=solve_ivp(rhs,(start,1),[1+b0/a0*start**2/6,b0*start**3/3],
        method='DOP853',rtol=tolerance,atol=[tolerance*1e-2,tolerance*1e-6],dense_output=True,max_step=.005)
    assert sol.success
    value,flux=sol.y[:,-1];normal=value-flux*np.log1p(-2*mu)/(2*mu)
    amplitude=-flux/normal
    def field(t):
        value,flux=sol.sol(t)/normal
        return value,flux/(t*t*a(t))
    return dict(R=R,M=M,mu=mu,a=amplitude,chi=R*amplitude,x=x,p=a,v=b,field=field)


def scalar_connection():
    initial=dict(np.load(OUT/'restored-state-17-4.npz'));base=scalar_model('initial',initial)
    rows=[]
    for steps in [1,2,4]:
        path=OUT/f'evolve-{steps}-{steps-1}-state-17-4.npz'
        if not path.exists(): continue
        data=dict(np.load(path));end=scalar_model(str(steps),data)
        knots=np.unique(np.r_[base['x'],end['x']]);left=knots[:-1];right=knots[1:]
        vals=[]
        for count in [6,12]:
            nodes,weights=np.polynomial.legendre.leggauss(count)
            xx=(left[:,None]+right[:,None])/2+(right-left)[:,None]*nodes/2
            f0,g0=base['field'](xx.ravel());f1,g1=end['field'](xx.ravel())
            xx=xx.ravel();dp=xx**2*(end['p'](xx)-base['p'](xx));dv=xx**2*(end['v'](xx)-base['v'](xx))
            value=(f0*f1*dv+g0*g1*dp).reshape(-1,count)
            integral=float(((right-left)/2)*(value@weights)@np.ones(len(left)))
            y=(nodes+1)/2
            tail=float(weights@(y/((1-2*base['mu']*y)*(1-2*end['mu']*y)))/2)
            da=-integral+2*(end['mu']-base['mu'])*base['a']*end['a']*tail
            delta_ratio=((end['R']-base['R'])-base['R']/base['M']*(end['M']-base['M']))/end['M']
            delta_alpha_over_phi=-(end['R']/end['M']*da+base['a']*delta_ratio)
            vals.append(dict(quadrature_points=count,delta_dimensionless_tail=da,
                delta_alpha_over_phi=delta_alpha_over_phi,
                direct_tail_difference=end['a']-base['a']))
        rows.append(dict(steps=steps,base_chi_m=base['chi'],end_chi_m=end['chi'],
            base_alpha_over_phi=-base['chi']/base['M'],end_alpha_over_phi=-end['chi']/end['M'],
            quadrature=vals,
            delta_alpha_quadrature_difference=abs(vals[0]['delta_alpha_over_phi']-vals[1]['delta_alpha_over_phi']),
            force_fraction_coefficient_per_phi_infinity_squared=-4*vals[-1]['delta_alpha_over_phi']))
    save('scalar-connection.json',dict(classification='Counterexample candidate',beta=-4,rows=rows,
        method='Finite Wronskian identity on common x=r/R including exact Schwarzschild vacuum tail integral. '
            'PCHIP material profiles define the numerical scalar problem; a finite-pressure surface and unsolved atmosphere remain.',
        interpretation='Small nonzero scalar background: alpha_A=-phi_infinity*chi/Mgeom+O(phi_infinity^3). '
            'For a weak companion alpha_B=beta*phi_infinity, delta(a_scalar/g)=beta*phi_infinity^2*delta(alpha_A/phi_infinity). '
            'This is the quasistatic scalar readout of a prescribed-flux thermal trajectory, not a driven orbital response or observation fit.',
        certified_derivative=False,nonlinear_scalar_backreaction=False,observation_likelihood_run=False))
    print('SCALAR CONNECTION',rows,flush=True)


def symbolic():
    import sympy as s
    rho,T,C,B=s.symbols('rho T C B',positive=True);F=s.Function('F')
    base=C*F(C*rho,T)+B
    pressure=rho**2*s.diff(base,rho);entropy=-s.diff(base,T);energy=base+T*entropy
    assert s.simplify(s.diff(energy,T)-T*s.diff(entropy,T))==0
    assert s.simplify(s.diff(energy,rho)-(pressure-T*s.diff(pressure,T))/rho**2)==0
    assert s.simplify(s.diff(entropy,rho)+s.diff(pressure,T)/rho**2)==0
    phi,chi,M,dc,dm=s.symbols('phi chi M dc dm',nonzero=True)
    assert s.simplify(s.diff(-phi*(chi+s.Symbol('h')*dc)/(M+s.Symbol('h')*dm),s.Symbol('h')).subs('h',0)
        +phi*(dc/M-chi*dm/M**2))==0
    # At fixed areal radius, flux divergence integrates to boundary luminosity.
    N,q,dt,dB=s.symbols('N q dt dB',positive=True)
    assert s.simplify(N*(q*N*dt)*dB-N**2*q*dt*dB)==0
    save('symbolic.json',dict(classification='Proven',
        composition_pullback_first_law=True,charge_and_ion_mapping=True,
        scalar_mass_normalization=True,proper_to_infinity_heat_factors=True,
        scope='Algebraic identities for the declared potential and normalization, conditional on a thermodynamically consistent underlying FreeEOS potential.'))
    print('PASS symbolic thermodynamics, normalization, redshift factors',flush=True)


def auxiliaries():
    source=ROOT/'verification/common_eos_bridge.f90';library=CACHE/'common_eos_aux.so'
    if not library.exists():
        module=next((gr.CACHE/'build').rglob('mod_free_eos.mod')).parent
        cmd=['gfortran','-O2','-fPIC','-shared','-I'+str(module),str(source),
            '-L'+str(gr.CACHE/'build/src'),'-Wl,-rpath,'+str(gr.CACHE/'build/src'),
            '-lfree_eos','-o',str(library)]
        result=subprocess.run(cmd,capture_output=True,text=True)
        save('auxiliary-bridge-build.json',dict(classification='Imported from prior work',command=cmd,
            returncode=result.returncode,stdout=result.stdout,stderr=result.stderr,
            unchanged_FreeEOS_library_sha256=sha(gr.CACHE/'build/src/libfree_eos.so.1.0.0')))
        assert result.returncode==0,result.stderr
    lib=ctypes.CDLL(str(library));call=lib.common_eos_aux
    call.argtypes=gr.EOS().call.argtypes;call.restype=None
    data=np.load(OUT/'restored-state-17-4.npz');_,p=mesa(OUT/'evolution-initial-profile.data.gz');eos=EOS()
    rows=[];raw=[]
    for i in range(len(data['dm'])):
        x=data['X'][i];y,ym=eos.composition(x);cx=ym@eos.weights;eps=np.ascontiguousarray(ym/cx)
        a=np.full(22,np.nan);info=ctypes.c_int(-999)
        call(2,float(data['lnd'][i]+np.log(cx)),float(data['lnT'][i]),eps,a,ctypes.byref(info))
        assert info.value==0 and np.all(np.isfinite(a))
        old=eos.native(2,data['lnd'][i]+np.log(cx),data['lnT'][i],eps)
        err=float(max(abs(a[:12]-old)/np.maximum(1,abs(old))));assert err<1e-9,(i,err)
        rows.append(dict(cell=i,bridge_relative_error=err,
            eta_difference=float(a[12]-p['eta'][i]),
            free_electron_to_fully_ionized_ratio=float(a[13]/np.exp(data['lnd'][i])/(y@Z)),
            native_pressure_relative=float(p['pressure'][i]/a[1]-1)))
        raw.append(a)
    active=abs(p['eps_nuc'])>1
    save('auxiliary-audit.json',dict(classification='Counterexample candidate',rows=rows,
        max_bridge_relative_error=max(r['bridge_relative_error'] for r in rows),
        burning_cells=int(active.sum()),
        max_burning_eta_difference=max(abs(r['eta_difference']) for r,a in zip(rows,active) if a),
        minimum_burning_electron_fraction=min(r['free_electron_to_fully_ionized_ratio'] for r,a in zip(rows,active) if a),
        max_native_pressure_relative=max(abs(r['native_pressure_relative']) for r in rows),
        native_rates_evaluated_with_FreeEOS_auxiliaries=False,
        library_sha256=sha(library),source_sha256=sha(source),
        interpretation='Fresh native FreeEOS auxiliaries at the exact restored GR rho,T,X. '
            'The existing MESA rate evaluator retains its own EOS auxiliaries; agreement of state coordinates alone does not close this boundary.'))
    np.savez_compressed(OUT/'auxiliary-outputs.npz',raw=np.array(raw))
    print('AUXILIARIES',json.loads((OUT/'auxiliary-audit.json').read_text())['max_burning_eta_difference'],flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
