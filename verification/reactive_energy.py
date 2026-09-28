"""Request24: composition-dependent energy references, not a reaction-rate solver.

Proven: reference changes leave total energy invariant with their composition term.
Counterexample candidate: prescribed finite burns using the frozen thermal EOS.
"""
from pathlib import Path
from decimal import Decimal as D
import json, re, shutil, sys
import gzip
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from thermal_restart import sha
from thermal_wd import mesa
import gr_mass as gr
import baryon_entropy as be

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/reactive-energy24'
OLD=ROOT/'outputs/baryon-entropy23'
MESA=Path('/home/lpaiu/work/thermal-restart19/mesa-r7624')
C2=(gr.C*100)**2
# Frozen MESA r7624 constants; Qconv converts MeV per reaction times mol/g/s.
QCONV=float(D('1.602176487e-6')*D('6.02214179e23'))
DOCS=['model-definition','observable-targets','adiabatic-limit',
      'nonadiabatic-regime','failure-ledger-dynamic-chi']


def save(name,obj):
    (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists()
    OUT.mkdir(parents=True)
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='9ab3723',
        previous_manifest_sha256=sha(OLD/'manifest.json'),
        primary='Keep Request23 neutral-atom rest mass and EOS energy unchanged. Audit the exact distributed network and mass conventions, then connect the nuclear-Q reference by an explicit composition-dependent internal-energy offset.',
        tests=['Decimal reaction stoichiometry and source Q reconstruction',
               'CNO cycle versus net 4H to He closure and electron counting',
               'Direct FreeEOS fixed-baryon-density composition derivatives',
               'Finite prescribed CNO burns in both energy references',
               'Negative controls omitting the offset and double-counting nuclear energy',
               'Symbolic first law, chemical term and GR luminosity identity'],
        tolerances=dict(reference_energy_relative=1e-10,finite_burn_energy_relative=1e-9,
            temperature_reference_relative=1e-10,derivative_refinement_relative=1e-4),
        boundary='Prescribed reaction extents have no time or rate. Source standard Q values are not runtime weak-rate Q evaluations. No complete network, transport, neutrino evolution, new TOV fit or observational inference is claimed.'))
    names=['chem/public/chem_lib.f90','chem/public/chem_def.f90',
        'chem/private/chem_isos_io.f90','const/public/const_def.f90',
        'rates/private/rates_initialize.f90','net/private/net_eval.f90',
        'net/private/net_derivs.f90','net/private/net_derivs_support.f90',
        'data/rates_data/reactions.list','data/rates_data/weak_info.list',
        'rates/public/rates_def.f90','data/chem_data/isotopes.data',
        'data/net_data/nets/basic.net','data/net_data/nets/add_cno_extras',
        'data/net_data/nets/add_hot_cno','star/private/eps_grav.f90']
    for name in names:
        dest=OUT/'sources'/name;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(MESA/name,dest)
    shutil.copy2(ROOT/'outputs/thermal-closure22/inputs/cno_extras.net',OUT/'sources/cno_extras.net')
    from thermal_robustness import gzcopy
    gzcopy(MESA/'data/rates_data/jina_reaclib_results_20130213default2',OUT/'sources/reaclib.gz')
    save('source-bindings.json',dict(classification='Imported from prior work',release=7624,
        sha256={p.relative_to(OUT/'sources').as_posix():sha(p) for p in sorted((OUT/'sources').rglob('*')) if p.is_file()},
        references=[dict(title='MESA r22.11.1: composition term in eps_grav',
            url='https://docs.mesastar.org/en/22.11.1/reference/controls.html#include-composition-in-eps-grav',
            use='Independent explanation of why composition derivatives matter; this later version is not evidence of r7624 runtime behavior.'),
            dict(title='MESA IV (2018)',url='https://arxiv.org/abs/1710.08424',
                 use='Primary energy-conservation context; numerical claims here use the archived r7624 sources.')]))


def isotopes():
    lines=(OUT/'sources/data/chem_data/isotopes.data').read_text().splitlines()[1:]
    result={}
    for line in lines[::4]:
        p=line.split();w,z,n,ex=D(p[1]),int(p[2]),int(p[3]),D(p[5])
        # chem_isos_io sets binding=0 for A<=1, hence uses del_Mp for H1.
        qex=D('7.288969') if p[0]=='h1' else D('8.071323') if p[0]=='neut' else ex
        result[p[0]]=dict(w=w,a=z+n,z=z,ex=ex,qex=qex)
    return result


def network():
    active=set();species=set()
    def visit(path):
        text='\n'.join(line.split('!')[0] for line in path.read_text().splitlines())
        pattern=r"include\s+'([^']+)'|(add_isos|add_reactions|remove_reactions|remove_reaction)\s*\(([^)]*)\)"
        for match in re.finditer(pattern,text):
            inc,op,body=match.groups()
            if inc: visit(OUT/'sources/data/net_data/nets'/inc)
            else:
                tokens=body.replace(',',' ').split()
                if op=='add_isos': species.update(tokens)
                elif op=='add_reactions': active.update(tokens)
                else: active.difference_update(tokens)
    visit(OUT/'sources/cno_extras.net')
    assert len(active)==72 and len(species)==22,(len(active),len(species))
    expected=json.loads((ROOT/'outputs/thermal-restart19/network-diagnosis.json').read_text())
    assert species==set(expected['profile_isotopes'])
    return active,species


def reactions():
    active,species=network();iso=isotopes();rows=[]
    energy_unit=D('29979245800')**2/(D('1.602176487e-6')*D('6.02214179e23'))
    for line in (OUT/'sources/data/rates_data/reactions.list').read_text().splitlines():
        name=line[:35].strip()
        if name not in active: continue
        nu={}
        for sign,field in [(-1,line[35:61]),(1,line[64:90])]:
            for count,key in re.findall(r'(\d+)\s+([a-z]+\d*)',field):
                nu[key]=nu.get(key,0)+sign*int(count)
        nu={k:v for k,v in nu.items() if v}
        qfield=line[105:123].strip().lower().replace('d','e')
        nfield=line[123:139].strip().lower().replace('d','e')
        baryon=sum(v*iso[k]['a'] for k,v in nu.items())
        complete=bool(nu) and baryon==0 and set(nu)<=species
        qstd=D(qfield) if qfield else -sum(v*iso[k]['qex'] for k,v in nu.items())
        qw=-sum(v*iso[k]['w'] for k,v in nu.items())*energy_unit if complete else None
        rows.append(dict(name=name,source='reactions.list',nu=nu,stoichiometry_complete=bool(complete),
            baryon_residual=baryon,explicit_Q_override=bool(qfield),
            Q_standard_MeV=str(qstd),Q_neutrino_list_MeV=str(D(nfield or '0')),
            Q_from_W_MeV=None if qw is None else str(qw),
            W_minus_standard_eV=None if qw is None else float((qw-qstd)*D('1e6'))))
    missing=active-{r['name'] for r in rows}
    weak={}
    for line in (OUT/'sources/data/rates_data/weak_info.list').read_text().splitlines():
        p=line.split()
        if len(p)>=4 and p[0] in iso and p[1] in iso:
            weak[(p[0],p[1])]=D(p[3].lower().replace('d','e'))
    with gzip.open(OUT/'sources/reaclib.gz','rt') as f: lines=f.read().splitlines()
    assert len(lines)%4==0
    for name in sorted(missing):
        _,left,kind,right=name.split('_')
        assert kind in ['wk','pg','ag','ga'],name
        inputs=[left]+(['h1'] if kind=='pg' else ['he4'] if kind=='ag' else [])
        outputs=[right]+(['he4'] if kind=='ga' else [])
        matching=[]
        for j in range(0,len(lines),4):
            chapter=int(lines[j]);parts=lines[j+1][5:35].split()
            parts=[{'p':'h1','n':'neut'}.get(k,k) for k in parts]
            if chapter!={('wk'):1,'pg':4,'ag':4,'ga':2}[kind]: continue
            if sorted(parts[:len(inputs)])==sorted(inputs) and sorted(parts[len(inputs):])==sorted(outputs):
                matching.append(j+1)
        assert matching,('no actual reaclib stoichiometry',name)
        nu={}
        for sign,keys in [(-1,inputs),(1,outputs)]:
            for key in keys: nu[key]=nu.get(key,0)+sign
        assert sum(v*iso[k]['a'] for k,v in nu.items())==0
        qstd=-sum(v*iso[k]['qex'] for k,v in nu.items())
        qw=-sum(v*iso[k]['w'] for k,v in nu.items())*energy_unit
        qnu=weak[(left,right)] if kind=='wk' else D(0)
        rows.append(dict(name=name,source='reaclib -> set_reaction_info -> get_Qtotal',
            reaclib_first_lines_1based=matching,nu=nu,stoichiometry_complete=True,
            baryon_residual=0,explicit_Q_override=False,Q_standard_MeV=str(qstd),
            Q_neutrino_list_MeV=str(qnu),Q_from_W_MeV=str(qw),
            W_minus_standard_eV=float((qw-qstd)*D('1e6'))))
    assert len(rows)==72 and len({r['name'] for r in rows})==72
    return sorted(rows,key=lambda r:r['name'])


def source_audit():
    rows=reactions();complete=[r for r in rows if r['stoichiometry_complete']]
    iso=isotopes();unit=D('29979245800')**2/(D('1.602176487e-6')*D('6.02214179e23'))
    reference_residuals=[]
    for row in complete:
        assert not row['explicit_Q_override']
        shift=-sum(v*((iso[k]['w']-iso[k]['a'])*unit-iso[k]['qex']) for k,v in row['nu'].items())
        reference_residuals.append(abs(D(row['Q_from_W_MeV'])-D(row['Q_standard_MeV'])-shift))
    assert max(reference_residuals)<D('1e-20')
    # Source-level standard-Q audit does not substitute for runtime weak Qs.
    ev=(OUT/'sources/net/private/net_eval.f90').read_text()
    assert 'actual_Qs(i) = n% Q(weak_id)' in ev
    ds=(OUT/'sources/net/private/net_derivs_support.f90').read_text()
    assert 'Q = reaction_Q - Qneu' in ds
    assert 'n% eps_neu_total = n% eps_neu_total + Qneu*rvs(i_rate)' in ds
    save('reaction-audit.json',dict(classification='Proven',network_species=22,network_entries=72,
        complete_stoichiometry_entries=len(complete),
        maximum_Decimal_reference_conversion_residual_MeV=str(max(reference_residuals)),
        hydrogen_mass_excess_vs_binding_override_eV=float((iso['h1']['ex']-iso['h1']['qex'])*D('1e6')),
        maximum_abs_W_minus_standard_eV=max(abs(r['W_minus_standard_eV']) for r in complete),
        unclosed_entries=[r['name'] for r in rows if not r['stoichiometry_complete']],
        interpretation='Auxiliary/special entries lacking a closed 22-isotope stoichiometry remain unaudited. Runtime weak Q and neutrino tables may override these standard values. Atomic Q already includes electron rest bookkeeping; do not add positron annihilation again.',rows=rows))
    print('Source entries',len(rows),'closed',len(complete),'max Q difference eV',max(abs(r['W_minus_standard_eV']) for r in complete))


def cycle_check():
    rows={r['name']:r for r in reactions()}
    names=['r_c12_pg_n13','r_n13_wk_c13','r_c13_pg_n14',
           'r_n14_pg_o15','r_o15_wk_n15','r_n15_pa_c12']
    nu={}
    for name in names:
        r=rows[name];assert r['stoichiometry_complete'] and not r['explicit_Q_override']
        for k,v in r['nu'].items(): nu[k]=nu.get(k,0)+v
    nu={k:v for k,v in nu.items() if v};assert nu=={'h1':-4,'he4':1},nu
    q=sum(D(rows[n]['Q_standard_MeV']) for n in names)
    qw=sum(D(rows[n]['Q_from_W_MeV']) for n in names)
    qnu=sum(D(rows[n]['Q_neutrino_list_MeV']) for n in names)
    iso=isotopes();assert q==4*iso['h1']['qex']-iso['he4']['qex']
    result=dict(classification='Proven',reactions=names,net_stoichiometry=nu,
        total_standard_MeV=float(q),total_W_MeV=float(qw),list_neutrino_MeV=float(qnu),
        standard_deposited_MeV=float(q-qnu),W_deposited_MeV=float(qw-qnu),
        difference_eV=float((qw-q)*D('1e6')),
        neutrino_reference='Prescribed sum of the two list mean energies, not re-evaluated neutrino spectra or rates.')
    save('cycle-control.json',result);print(json.dumps(result,indent=2))


class Energy:
    def __init__(self):
        self.iso=isotopes()
        self.keys=[k for k in json.loads((ROOT/'outputs/thermal-restart19/network-diagnosis.json').read_text())['profile_isotopes'] if not k.startswith('f')]
        self.A=np.array([self.iso[k]['a'] for k in self.keys])
        self.W=np.array([float(self.iso[k]['w']) for k in self.keys])
        self.EX=np.array([float(self.iso[k]['qex']) for k in self.keys])
        self.elements=np.array([gr.ELEMENTS.index(re.sub(r'\d','',k)) for k in self.keys])
        self.eos=gr.EOS()
        # Stable excess-energy form avoids subtracting two ~c^2 terms.
        self.g_per_Y=(self.W-self.A)*C2-self.EX*QCONV

    def state(self,X,rhoB,lt):
        X=np.asarray(X);assert len(X)==len(self.A) and np.all(X>=0) and abs(sum(X)-1)<1e-12
        y=X/self.A;cx=np.dot(y,self.W)
        eps=np.bincount(self.elements,weights=y,minlength=20)/cx
        a=self.eos(2,np.log(rhoB*cx),lt,eps)
        u=cx*a[2];s=cx*a[3];g=np.dot(y,self.g_per_Y)
        return dict(cx=cx,u=u,s=s,g=g,uQ=u+g,p=a[1],cv=cx*a[10]/np.exp(lt))


def material_states():
    e=Energy();m=be.Material();_,d=mesa(ROOT/'outputs/thermal-closure22/selected.data.gz')
    X=np.array([d[k] for k in e.keys]).T;X/=X.sum(axis=1)[:,None]
    peak=int(np.argmax(d['eps_nuc']))
    # Reproduce the scaled Request23 branch; keep its frozen data untouched.
    solver=be.Solver(points=33,subdivision=4)
    row=json.loads((OLD/'scaled-33-4.json').read_text())
    params=[row['parameters'][0],row['parameters'][1],0.]
    error,inner,outer=solver.branches(params,row['baryon_scale'],record=True)
    assert max(abs(error))<1e-8
    # Cell-face interpolation in each branch's material coordinate. Samples
    # are midpoint diagnostics, not a new hydrostatic solution or mesh limit.
    records=[]
    for i in sorted(set([0,len(X)//2,len(X)-1,peak])):
        if i>=m.split:
            q=(m.inner[i]+m.inner[i+1])/2
            lp=np.interp(q,inner[:,0],inner[:,3])
        else:
            w=(m.outer[i]+m.outer[i+1])/2
            lp=np.interp(w,outer[:,0],outer[:,3])
        a,lt,_=be.invert(e.eos,lp,solver.ref[i,3],m.eps[i],m.lt[i])
        rhoB=a[0]/m.cx[i]
        records.append(dict(zone=int(d['zone'][i]),index=i,X=X[i].tolist(),rhoB=rhoB,logT=lt,
            source_nuclear_peak=i==peak,origin='Scaled Request23 material cell midpoint, entropy-inverted FreeEOS'))
    return records


def eos_controls():
    e=Energy();records=material_states();rows=[]
    direction=np.zeros(len(e.keys));direction[e.keys.index('h1')]=-1;direction[e.keys.index('he4')]=1
    for rec in records:
        X=np.array(rec['X']);rho=rec['rhoB'];lt=rec['logT'];a=e.state(X,rho,lt)
        # H burning direction requires non-negligible H on both finite-difference sides.
        if min(X[e.keys.index('h1')],X[e.keys.index('he4')])<1e-8: continue
        step=min(1e-5,.01*min(X[e.keys.index('h1')],X[e.keys.index('he4')]))
        derivatives=[]
        for h in [step,step/2,step/4]:
            p=e.state(X+h*direction,rho,lt);n=e.state(X-h*direction,rho,lt)
            du=(p['u']-n['u'])/(2*h);tds=np.exp(lt)*(p['s']-n['s'])/(2*h)
            derivatives.append([du,tds,du-tds])
        last=np.array(derivatives[-1]);prev=np.array(derivatives[-2])
        diff=float(np.max(abs(last-prev)/np.maximum(abs(last),1.)))
        assert diff<1e-4,(rec['zone'],diff)
        h=1e-4;sp=e.state(X,rho,lt+h);sn=e.state(X,rho,lt-h)
        cvfd=(sp['u']-sn['u'])/(2*h*np.exp(lt))
        firstlaw=(sp['u']-sn['u'])/(np.exp(lt)*(sp['s']-sn['s']))-1
        assert abs(cvfd/a['cv']-1)<1e-5 and abs(firstlaw)<1e-5
        rows.append(dict(zone=rec['zone'],rhoB=rho,T_K=float(np.exp(lt)),step=step,
            fixed_rhoB_T_derivatives_per_H_fraction=derivatives,
            last_refinement_relative=diff,fixed_X_firstlaw_relative=float(firstlaw),
            heat_capacity_FD_relative=float(cvfd/a['cv']-1),
            thermal_chemical_term_over_Tds=float(last[2]/last[1])))
    assert len(rows)>=2
    save('material-samples.json',dict(classification='Counterexample candidate',records=records))
    save('eos-composition-controls.json',dict(classification='Counterexample candidate',results=rows,
        interpretation='du_B - T ds_B at fixed rho_B,T is the chemical composition term. No reaction rate, physical continuum bound or actual MESA absolute-entropy match is inferred.'))
    print('EOS composition controls',len(rows),'max refinement',max(r['last_refinement_relative'] for r in rows))


def burn_controls():
    e=Energy();rec=next(r for r in json.loads((OUT/'material-samples.json').read_text())['records'] if r['source_nuclear_peak'])
    X0=np.array(rec['X']);rho=rec['rhoB'];lt0=rec['logT'];a0=e.state(X0,rho,lt0)
    cyc=json.loads((OUT/'cycle-control.json').read_text());qnu=cyc['list_neutrino_MeV']
    direction=np.zeros(len(e.keys));direction[e.keys.index('h1')]=-1;direction[e.keys.index('he4')]=1
    amount=min(1e-4,.001*X0[e.keys.index('h1')]);rows=[]
    for burn in np.linspace(amount/32,amount,32):
        dx=burn*direction;dy=dx/e.A;X=X0+dx;extent=burn/4
        restchange=float(np.dot(dy,e.W-e.A)*C2)
        nu_loss=extent*qnu*QCONV;heatW=-restchange-nu_loss
        heatQ=-float(np.dot(dy,e.EX)*QCONV)-nu_loss
        def solve(field,target):
            return brentq(lambda lt:e.state(X,rho,lt)[field]-target,lt0-.7,lt0+.7,xtol=2e-14)
        ltW=solve('u',a0['u']+heatW);a=e.state(X,rho,ltW)
        ltQ=solve('uQ',a0['uQ']+heatQ)
        residual=(a['u']-a0['u']+restchange+nu_loss)/heatW
        refdiff=np.expm1(ltQ-ltW)
        assert abs(residual)<1e-9 and abs(refdiff)<1e-10
        # Deliberately wrong: old EOS u plus unconverted standard nuclear heat.
        ltbad=solve('u',a0['u']+heatQ);bad=e.state(X,rho,ltbad)
        badres=(bad['u']-a0['u']+restchange+nu_loss)/heatW
        # If total energy includes rest loss, adding the same heat as a source
        # again deposits twice. This is a separate explicit negative control.
        ltdouble=solve('u',a0['u']+2*heatW);double=e.state(X,rho,ltdouble)
        doubleres=(double['u']-a0['u']+restchange+nu_loss)/heatW
        rows.append(dict(H_fraction_burned=float(burn),extent_mol_g=float(extent),
            T_K=float(np.exp(ltW)),heat_W_erg_g=heatW,heat_Q_erg_g=heatQ,
            neutrino_erg_g=nu_loss,energy_relative=float(residual),temperature_reference_relative=float(refdiff),
            omitted_offset_energy_relative=float(badres),double_count_energy_relative=float(doubleres),
            entropy_change_erg_g_K=float(a['s']-a0['s'])))
    assert max(abs(r['omitted_offset_energy_relative']) for r in rows)>1e-7
    assert min(abs(r['double_count_energy_relative']) for r in rows)>.99
    # Finite reference shift of arbitrary magnitude, with total energy invariant.
    r=rows[-1];X=X0+r['H_fraction_burned']*direction;a=e.state(X,rho,np.log(r['T_K']))
    reference_change=(np.dot(X/e.A,e.EX)*QCONV+a['uQ'])-(np.dot(X/e.A,e.W-e.A)*C2+a['u'])
    assert abs(reference_change)<1e-10*max(abs(a['u']),1)
    save('finite-burn-controls.json',dict(classification='Counterexample candidate',zone=rec['zone'],
        rhoB=rho,T0_K=float(np.exp(lt0)),prescribed_neutrino_MeV=qnu,
        time_or_rate_assigned=False,hydrostatic_or_transport_evolution=False,
        reference_total_energy_difference_erg_g=float(reference_change),results=rows))
    print('Finite burn final',json.dumps(rows[-1],indent=2))


def differential_control():
    e=Energy();rec=next(r for r in json.loads((OUT/'material-samples.json').read_text())['records'] if r['source_nuclear_peak'])
    burn=json.loads((OUT/'finite-burn-controls.json').read_text())['results'][-1]
    X0=np.array(rec['X']);rho=rec['rhoB'];lt0=rec['logT'];amount=burn['H_fraction_burned']
    direction=np.zeros(len(e.keys));direction[e.keys.index('h1')]=-1;direction[e.keys.index('he4')]=1
    q=burn['heat_W_erg_g']/amount;rows=[]
    for h,chemical in [(1e-5,True),(5e-6,True),(5e-6,False)]:
        def rhs(z,y):
            X=X0+z*amount*direction;lt=float(y[0]);T=np.exp(lt)
            a=e.state(X,rho,lt);p=e.state(X+h*direction,rho,lt);n=e.state(X-h*direction,rho,lt)
            uX=(p['u']-n['u'])/(2*h);TsX=T*(p['s']-n['s'])/(2*h)
            term=uX if chemical else TsX
            return [amount*(q-term)/(T*a['cv'])]
        sol=solve_ivp(rhs,(0.,1.),[lt0],method='DOP853',rtol=1e-11,atol=1e-12,max_step=.1)
        assert sol.success
        T=float(np.exp(sol.y[0,-1]));a=e.state(X0+amount*direction,rho,np.log(T));a0=e.state(X0,rho,lt0)
        error=T/burn['T_K']-1;energy=(a['u']-a0['u'])/burn['heat_W_erg_g']-1
        if chemical: assert abs(error)<1e-7 and abs(energy)<1e-7
        else: assert abs(energy)>1e-4
        rows.append(dict(composition_difference_step=h,chemical_term_included=chemical,
            T_K=T,temperature_vs_finite_energy_relative=error,energy_relative=float(energy),calls=sol.nfev))
    assert abs(rows[1]['T_K']/rows[0]['T_K']-1)<1e-7
    save('differential-control.json',dict(classification='Counterexample candidate',results=rows,
        interpretation='Independent DOP853 integration over prescribed reaction extent, using EOS heat capacity and numerical composition partials. This is a constant-density reactor control, not time evolution or a certified derivative bound.'))
    print('Differential control',json.dumps(rows,indent=2))


def symbolic():
    import sympy as s
    u,r,g,nu,qext,p,v=s.symbols('udot rdot gdot nu qext p vdot')
    # Internal equation in W reference: udot + p vdot = -rdot - nu + qext.
    total=s.expand((u+p*v+r+nu-qext)-(u+p*v-(-r-nu+qext)))
    shifted=s.expand(((u+g)+p*v+(r-g)+nu-qext)-(u+p*v+r+nu-qext))
    T,sdot,muY,nuth,src,red=s.symbols('T sdot muY nuth src red')
    # Correct GR luminosity source with composition; q_nuc excludes reaction nu.
    luminosity=red*(src-nuth-T*sdot-muY)
    assert total==0 and shifted==0
    assert s.simplify(luminosity.subs(muY,u+p*v-T*sdot)-red*(src-nuth-u-p*v))==0
    assert s.diff(luminosity,muY)==-red
    save('symbolic.json',dict(classification='Proven',total_energy_identity=True,
        composition_reference_invariance=True,chemical_term_in_GR_luminosity=True,
        equations=['e_B = C_X c^2 + u_B',
            'du_B/dtau + P d(1/rho_B)/dtau = -c^2 dC_X/dtau - q_nu + q_ext',
            'u_Q = u_W + g(X); e0_Q = e0_W - g(X); q_Q = q_W + dg/dtau',
            'du_B + P dv_B = T ds_B + sum(mu_thermal_i dY_i)',
            'dL_inf/dB = exp(2nu) [q_nuc - q_thermal_nu - T ds_B/dtau - sum(mu_thermal_i dY_i/dtau)]'],
        assumptions='Comoving conserved baryon mass, electrically neutral equilibrium EOS, transparent emitted neutrinos, quasistatic GR luminosity balance as in Request22. No diffusive species flux; otherwise include chemical-energy flux. Rates must use the same reference as u.'))
    print('PASS symbolic total-energy, reference-shift and chemical-luminosity identities')


def maintain():
    prior=json.loads((OLD/'manifest.json').read_text())['sha256'];bindings={}
    dest=OUT/'request23-notes';dest.mkdir(exist_ok=False)
    for rel in ['docs/'+k+'.md' for k in DOCS]+['paper/revision-manifest.json']:
        assert sha(ROOT/rel)==prior[rel],rel
        snap=dest/Path(rel).name;shutil.copy2(ROOT/rel,snap)
        bindings[rel]=dict(snapshot=snap.relative_to(ROOT).as_posix(),sha256=prior[rel],historical_manifest=rel.startswith('paper/'))
    save('historical-note-bindings.json',bindings)
    paragraphs=[
        '분류: Proven. Request24는 조성 의존 에너지 기준을 연결했다. e_B=C_X c²+u_W를 유지하면서 u_Q=u_W+g(X), e0_Q=e0_W−g(X)로 바꾸면 총 에너지는 그대로다. g는 배포 원자량과 MESA 표준 Q를 만드는 질량초과의 차이이다. 따라서 핵 가열항도 dg/dτ만큼 함께 바뀌어야 한다. 22종 네트워크의 72항목 중 닫힌 반응식 64개에서 이 변환을 검산했으며 나머지 8개 보조율은 독립 반응식 인증 대상에서 제외했다. 이 기준 변환은 Request23의 GR 질량을 다시 적합하거나 바꾸지 않는다.',
        '분류: Counterexample candidate. 지정된 CNO 순반응량의 직접 EOS 시험에서 기준 보정 누락은 열수지 6.21335 ppm, 화학적 조성 항을 뺀 엔트로피 적분은 0.481113% 오차를 냈다. 보정한 유한 에너지 역산과 독립 미분 적분은 일치한다. 이들은 에너지 장부 검증이며 실제 반응 시간, 궤도 완화시간이나 새로운 관측량은 아니다.',
        '분류: Proven. 조성이 변하면 du_B+P dv_B=T ds_B+Σ μ_i^th dY_i이다. 따라서 ds_B=0만으로 du_B+P dv_B=0을 결론낼 수 없다. 고정 조성의 Request23 엔트로피 재구성은 보존되지만, 반응이 켜진 후 동일한 엔트로피를 고정하는 것은 일반적인 열 진화 해가 아니다.',
        '분류: Proven. Request22의 준정적 GR 광도 식을 반응하는 물질로 확장할 때 dL∞/dB=e^(2ν)[q_nuc−q_thermalν−T ds_B/dτ−Σ μ_i^th dY_i/dτ]로 화학 항을 명시해야 한다. q_nuc는 같은 에너지 기준에서 반응 중성미자를 이미 뺀 가열률이다. 바리온과 함께 움직이는 물질, 종별 확산 유속 없음, 투명한 중성미자 가정이며 확산 시 화학 에너지 유속도 필요하다.\n\n분류: Conjectural. 다음 단계는 새 상태에서의 실제 반응률·중성미자·복사/전도/대류와 GR 열·조성 진화이다. 현재 시험은 반응량을 지정했으며 속도나 시간을 계산하지 않았다.',
        '분류: Proven. 조성 변화에 고정 조성 T ds 항등식을 그대로 적용하거나, 총 에너지에 정지질량 감소를 포함하면서 핵 가열을 다시 더하는 연결은 실패한다. 전자는 화학적 조성 항, 후자는 일관된 내부/총 에너지 식 선택이 필요하다. 배포 원자량과 표준 Q 표의 작은 차이도 조성 의존 기준 보정으로 명시했다.\n\n분류: Counterexample candidate. 19종 FreeEOS 조성으로 지정 CNO 반응량의 국소 검증을 통과했다. 22종 전체 EOS, 실제 weak Q/중성미자 재평가와 72항목 네트워크 시간 적분은 미완료이다. 모든 prior 실패, 원래 물질 보존 질량 불일치 및 별도 바리온 배율 후보 구분을 유지한다.'
    ]
    revision=json.loads((ROOT/'paper/revision-manifest.json').read_text())
    for stem,body in zip(DOCS,paragraphs):
        rel='docs/'+stem+'.md'
        with (ROOT/rel).open('ab') as f:
            f.write(('\n\n## Request 24 조성 변화의 반응 에너지 기준\n\n'+body+
                '\n\n세부 근거: [한글 보고서](../notes/REQUEST24_REACTIVE_ENERGY_KO.md).\n').encode())
        revision['sha256'][rel]=sha(ROOT/rel)
    revision['request24_supporting_note_update']=dict(evidence_manifest='outputs/reactive-energy24/manifest.json',
        historical_notes='outputs/reactive-energy24/historical-note-bindings.json',
        status='정지질량·핵 Q·내부에너지 기준 변환 및 조성 화학 항 검증; 실제 반응률·GR 열 진화 미완료',
        artifact_status='원고 PDF와 ZIP은 Request 12 동결본이다.')
    (ROOT/'paper/revision-manifest.json').write_text(json.dumps(revision,ensure_ascii=False,indent=2)+'\n')


def seal():
    import scipy,sympy
    sources=json.loads((OUT/'source-bindings.json').read_text())
    sources['sha256']={p.relative_to(OUT/'sources').as_posix():sha(p) for p in sorted((OUT/'sources').rglob('*')) if p.is_file()}
    save('source-bindings.json',sources)
    save('gates.json',dict(classification='Proven',standard_closed_reaction_references_reconciled=True,
        EOS_composition_chemical_term_checked=True,finite_burn_reference_invariance_passed=True,
        independent_differential_energy_control_passed=True,
        prior_GR_mass_or_material_results_changed=False,actual_reaction_rates_evaluated=False,
        full_22_isotope_reactive_EOS=False,GR_thermal_evolution_solved=False,
        global_derivative_error_certificate=False,complete_nonlinear_observational_inference=False,
        final_PDF_or_ZIP_updated=False,
        classification_detail='Theorem progress: energy-reference invariance and explicit chemical term in GR heat balance. Loophole progress: finite composition-change EOS controls; no new dynamic observable.'))
    save('provenance.json',dict(classification='Proven',before_task_checkpoint='9ab3723',
        previous_manifest_sha256=sha(OLD/'manifest.json'),interpreter=sys.executable,
        versions={m.__name__:m.__version__ for m in [np,scipy,sympy]},
        runtime='Unchanged FreeEOS bridge/runtime bound in outputs/gr-mass21/provenance.json; no MESA rebuild.',
        source_Q='Standard Q from mass_excess and A<=1 binding override, not the Reaclib fit Q field. Eleven auto-defined reactions have their stoichiometry checked against actual Reaclib records; weak_info.list supplies declared mean neutrino energy.',
        pilot_corrections=['Initial reactions.list-only parser could not resolve 11 auto-defined reactions; archived Reaclib and weak_info.list and checked them before source verdict.',
            'Initial negative-control temperature bracket +/-0.2 in lnT did not enclose double-heating roots. Expanded to +/-0.7 without changing reaction extents or acceptance thresholds. No failed physical verdict was replaced.',
            'Combined EOS recheck changed energy residuals by about 2e-15, so exact JSON float equality failed. Frozen results are retained; numeric reproduction uses 1e-12 times max(1,abs(reference)), with all original scientific acceptance gates rerun unchanged.'],
        sample_boundary='Four material midpoints from frozen scaled Request23; three with sufficient H/He support centred composition differences. Central hydrogen-depleted sample is not burned. Same 5735-cell template; not continuum convergence.'))
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    files += [ROOT/'verification/reactive_energy.py',ROOT/'notes/REQUEST24_REACTIVE_ENERGY_KO.md',ROOT/'.gitattributes']
    files += [ROOT/k for k in json.loads((OUT/'historical-note-bindings.json').read_text())]
    save('manifest.json',dict(classification='Proven',scope='Composition-dependent EOS/rest-energy/Q reference and prescribed-burn controls',
        sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)}))


def verify():
    stages=['validated-variational','remaining-levers15','nbody-readout16','nonzero-drive17',
        'thermal-wd18','thermal-restart19','thermal-robustness20','gr-mass21','thermal-closure22','baryon-entropy23','reactive-energy24']
    histories={f'outputs/{a}/manifest.json':f'outputs/{b}/historical-note-bindings.json' for a,b in zip(stages,stages[1:])}
    count=0
    for label in ['outputs/research-remediation/manifest.json',*histories,'paper/revision-manifest.json','outputs/reactive-energy24/manifest.json']:
        old=json.loads((ROOT/histories[label]).read_text()) if label in histories else {}
        for name,expected in json.loads((ROOT/label).read_text())['sha256'].items():
            path=ROOT/name
            if name in old:
                bind=old[name];path=ROOT/bind['snapshot'];assert expected==bind['sha256']
                if bind.get('historical_manifest'):
                    before=json.loads(path.read_text());after=json.loads((ROOT/name).read_text())
                    for k,value in before.items():
                        if k!='sha256': assert after[k]==value,k
                    for k,value in before['sha256'].items():
                        if k not in old: assert after['sha256'][k]==value,k
                else: assert (ROOT/name).read_bytes().startswith(path.read_bytes())
            assert sha(path)==expected,(label,name);count+=1
    gates=json.loads((OUT/'gates.json').read_text())
    for key in ['prior_GR_mass_or_material_results_changed','actual_reaction_rates_evaluated',
        'full_22_isotope_reactive_EOS','GR_thermal_evolution_solved','global_derivative_error_certificate',
        'complete_nonlinear_observational_inference','final_PDF_or_ZIP_updated']: assert not gates[key],key
    for key in ['runtime_sha256','module_sha256']:
        for path,digest in json.loads((ROOT/'outputs/gr-mass21/provenance.json').read_text())[key].items(): assert sha(path)==digest,path
    print('PASS:',count,'현재·역사 SHA 및 에너지/반응률/GR 진화 경계')


def recheck():
    import tempfile
    global OUT
    original=OUT
    worst=[0.]
    def compare(a,b,path):
        if isinstance(b,float):
            assert np.isfinite(a) and np.isfinite(b),path
            error=abs(a-b)/max(1.,abs(b));worst[0]=max(worst[0],error)
            assert error<=1e-12,(path,a,b,error)
        elif isinstance(b,dict):
            assert a.keys()==b.keys(),path
            for key in b: compare(a[key],b[key],path+'/'+key)
        elif isinstance(b,list):
            assert len(a)==len(b),path
            for i,(x,y) in enumerate(zip(a,b)): compare(x,y,path+'/'+str(i))
        else: assert a==b,(path,a,b)
    # Re-run into temporary output; never rewrite sealed numerical evidence.
    with tempfile.TemporaryDirectory(prefix='reactive-energy24-') as folder:
        OUT=Path(folder)/'out';shutil.copytree(original,OUT)
        try:
            source_audit();cycle_check();eos_controls();burn_controls();differential_control();symbolic()
            for name in ['reaction-audit.json','cycle-control.json','material-samples.json',
                'eos-composition-controls.json','finite-burn-controls.json','differential-control.json','symbolic.json']:
                compare(json.loads((OUT/name).read_text()),json.loads((original/name).read_text()),name)
        finally: OUT=original
    print('PASS seven numerical/symbolic/source outputs; maximum normalized reproduction difference',worst[0])


if __name__=='__main__': globals()[sys.argv[1]]()
