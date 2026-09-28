"""Audit physical partition prerequisites without replacing a native EOS library."""
import ctypes, gzip, io, json, math, re, shutil, sys
import numpy as np
import mpmath as mp
import sympy as sp
import eos_molecular_switch as switch

g=switch.g
DATA=g.OUT/'gr-molecular-partition-source'
OUT=g.OUT/'gr-molecular-partition-audit'
SOURCE=switch.CACHE/'source/src/mod_molecular_hydrogen.f90'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();switch.verify()
    shutil.copy2(SOURCE,OUT/SOURCE.name)
    files={'h2-2016.html':'https://arxiv.org/html/1607.04479v1',
        'cds-ReadMe.txt':'https://cdsarc.cds.unistra.fr/ftp/J/A+A/595/A130/ReadMe',
        'cds-index.html':'https://cdsarc.cds.unistra.fr/ftp/J/A+A/595/A130/',
        'h2equil.dat.gz':'https://cdsarc.cds.unistra.fr/ftp/J/A+A/595/A130/h2equil.dat.gz',
        'eleroyh2.dat':'http://kurucz.harvard.edu/molecules/h2/eleroyh2.dat'}
    for name in files:assert (DATA/name).stat().st_size>0
    raw=gzip.decompress((DATA/'h2equil.dat.gz').read_bytes())
    assert len(raw.splitlines())==20000 and len(raw)==3560000
    save('sources.json',dict(classification='Imported from prior work',retrieved_date='2026-09-11',
        files={name:dict(url=url,sha256=g.c.sha(DATA/name)) for name,url in files.items()},
        gzip_uncompressed_bytes=len(raw),gzip_crc_verified=True,
        cds='Popovas and Jorgensen 2016, A&A 595 A130, DOI 10.1051/0004-6361/201527209. Equilibrium H2, 1 K steps, 1--20000 K. Units/columns from the archived CDS ReadMe; the Eint/RT column is dimensionless by article equation 40 despite the catalogue unit label J.',
        levels='The 348 (v,J,term value) rows at the Kurucz URL cited by article footnote 4/Figure 1. Ground-electronic-state level list includes high levels; do not certify it as a complete bound spectrum. Term values are in cm^-1 relative to its (0,0) row. Use every listed row, with normalized equilibrium nuclear-spin weights (2J+1)/4 for even J and 3(2J+1)/4 for odd J, as in article equation 37 for the ground state.',
        transport='CDS HTTPS gzip completed with HTTP 200 then a successful 206 byte-range resume; gzip CRC checked. Kurucz data retrieved from the article original HTTP URL; HTTPS certificate mismatch was not bypassed. HTTP transport does not authenticate provenance; the archived SHA pins these received bytes only.',
        unavailable='UGA H2+ level URL linked from its mainh2p page timed out via HTTPS and HTTP. No H2+ spectral level data were obtained or substituted with HD+.',
        boundary='Published computed tables and spectral lists are numerical source data, not rigorous physical uncertainty enclosures, plasma occupation probabilities or a high-temperature complete EOS.'))
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='9de5955',
        discovery='The preliminary direct partition call at 1e6 K returned q_lnT+q_lnT_lnT < 0 for both molecules. The following audit is confirmation of that discovered failure, not a blind preregistered discovery.',
        comparison_T_K=[1000,3000,5000,8999,9000,10000,15000,18000,19000,20000,100000,1000000],
        table_comparison='Every integer temperature 1000--20000 K; no fit, normalization adjustment, rejection gate or extrapolated published table. Record native versus source differences and the source column identities.',
        native_replay_absolute_tolerance=1e-10,
        uniform_T_K=[999999,1000001],interval_digits=70,maximum_logT_derivative=10,
        printed_level_halfwidth_cm_inverse='0.0005',
        interval_scope='Outward interval calculation of the explicitly declared finite spectral sum and rounded-coefficient Taylor model. The level halfwidth is a printing-bin scenario, not an experimental/theoretical uncertainty bound. The polynomial interval excludes native compiler/libm evaluation error. These are not bounds on the implicit chemical-equilibrium EOS derivatives.',
        runtime=switch.read('manifest.json')['runtime'],
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
            g.ROOT/'verification/molecular_partition_data.py',OUT/'sources.json',OUT/SOURCE.name,
            g.OUT/'initial-state-17-4.npz',*[DATA/name for name in files]]}))


def coefficients():
    text=(OUT/SOURCE.name).read_text();result={}
    for name in ['qh2i90','qh2i90_hi','qh2plusi90','qh2plusi90_hi']:
        match=re.findall(r'\bdata '+name+r'/&(.*?)/',text,re.S);assert len(match)==1,name
        values=[x.strip() for x in match[0].replace('&','').replace('_fp_kind','').split(',')]
        result[name]=np.array([float(x) for x in values])
    return result


def native_function():
    eos=switch.EOS();f=eos.inventory_lib.__mod_molecular_hydrogen_MOD_molecular_hydrogen
    f.restype=None;f.argtypes=[ctypes.POINTER(ctypes.c_int)]*3+[ctypes.POINTER(ctypes.c_double)]*7
    def call(T):
        integers=[ctypes.c_int(i) for i in [0,3,2]];t=ctypes.c_double(np.log(float(T)))
        out=[ctypes.c_double() for _ in range(6)]
        f(*[ctypes.byref(i) for i in integers],ctypes.byref(t),*[ctypes.byref(i) for i in out])
        return np.array([v.value for v in out]).reshape(2,3)
    return call


def replay(T,coefs):
    base=min(float(1.000000000000001e5),max(900.,float(T)));x=np.log10(5040./base)
    d=np.log(T)-np.log(base);out=[]
    for name in ['qh2i90','qh2plusi90']:
        c=coefs[name+('_hi' if base>=9000 else '')]
        q=np.polynomial.polynomial.polyval(x,c)*np.log(10.)
        q1=-np.polynomial.polynomial.polyval(x,np.polynomial.polynomial.polyder(c))
        q2=np.polynomial.polynomial.polyval(x,np.polynomial.polynomial.polyder(c,2))/np.log(10.)
        out.append([q+d*q1+.5*d*d*q2,q1+d*q2,q2])
    return np.array(out)


def interval_record(value):
    return dict(display=str(value),binary_endpoints=[list(p) for p in value._mpi_])


def interval_audit(level_rows,coefs,plan):
    mp.iv.dps=plan['interval_digits'];iv=mp.iv
    def exact_float(x):
        n,d=float(x).as_integer_ratio();return iv.mpf(n)/d
    T=iv.mpf(plan['uniform_T_K']);T0=exact_float(1.000000000000001e5)
    x=iv.ln(5040/T0)/iv.ln(10);delta=iv.ln(T/T0);native_rows=[]
    def poly(c):
        value=iv.mpf(0)
        for a in reversed(c):value=value*x+a
        return value
    for name in ['qh2i90_hi','qh2plusi90_hi']:
        c=[exact_float(v) for v in coefs[name]]
        q1=-poly([i*c[i] for i in range(1,len(c))])
        q2=poly([i*(i-1)*c[i] for i in range(2,len(c))])/iv.ln(10)
        heat=q1+(delta+1)*q2
        assert heat.b<0
        native_rows.append(dict(name=name,internal_Cv_over_kB=interval_record(heat),strictly_negative=True))
    # D = d/d ln T; D^n exp(-x) = exp(-x) P_n(x), x = E/(k_B T).
    z=sp.symbols('x');polynomials=[sp.Integer(1)]
    for _ in range(plan['maximum_logT_derivative']):
        p=polynomials[-1];polynomials.append(sp.expand(z*(p-sp.diff(p,z))))
    c2=iv.mpf('6.62607015e-34')*299792458*100/iv.mpf('1.380649e-23')
    q=[iv.mpf(0) for _ in polynomials]
    for v,J,printed in level_rows:
        energy=iv.mpf(printed);width=iv.mpf(plan['printed_level_halfwidth_cm_inverse'])
        energy=iv.mpf([energy.a-width.b,energy.b+width.b]) if (v,J)!=(0,0) else iv.mpf(0)
        x=c2*energy/T;weight=iv.mpf((2*J+1)*(1 if J%2==0 else 3))/4
        boltzmann=weight*iv.exp(-x)
        for n,p in enumerate(polynomials):
            q[n]+=boltzmann*poly([int(v) for v in reversed(sp.Poly(p,z).all_coeffs())])
    assert q[0].a>0
    jets=[iv.ln(q[0])]
    for n in range(1,len(q)):
        jets.append((q[n]-sum(math.comb(n-1,k-1)*jets[k]*q[n-k] for k in range(1,n)))/q[0])
    heat=jets[1]+jets[2];assert heat.a>0
    return dict(classification='Proven',temperature_domain_K=plan['uniform_T_K'],
        rounded_coefficient_Taylor_model=native_rows,
        fixed_348_level_model=dict(logQ_derivatives=[interval_record(v) for v in jets],
            internal_Cv_over_kB=interval_record(heat),strictly_positive=True),
        mathematical_scope=plan['interval_scope'],mpmath_version=mp.__version__)


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in plan['runtime'].items():assert g.c.sha(path)==digest,path
    table=np.loadtxt(io.StringIO(gzip.decompress((DATA/'h2equil.dat.gz').read_bytes()).decode('ascii')))
    assert table.shape==(20000,9) and np.array_equal(table[:,0],np.arange(1,20001))
    assert np.all(np.isfinite(table)) and np.all(table[:,1]>0) and table[0,1]==.25
    rows=[]
    for line in (DATA/'eleroyh2.dat').read_text().splitlines():
        if line.strip():
            v,J,E=line.split();rows.append((int(v),int(J),E))
    assert len(rows)==348 and len({(v,J) for v,J,E in rows})==348
    assert rows[0]==(0,0,'0.000') and all(v>=0 and J>=0 and float(E)>=0 for v,J,E in rows)
    coefs=coefficients();native=native_function();temperatures=table[999:,0]
    values=np.array([native(T) for T in temperatures]);replayed=np.array([replay(T,coefs) for T in temperatures])
    error=float(np.max(abs(values-replayed)));assert error<plan['native_replay_absolute_tolerance'],error
    rho=np.load(g.OUT/'initial-state-17-4.npz');actual=np.exp(rho['lnT']);inside=np.flatnonzero(actual<=20000)
    selected=[]
    levels=np.array([float(E) for v,J,E in rows]);weights=np.array([(2*J+1)*(1 if J%2==0 else 3)/4 for v,J,E in rows])
    c2=6.62607015e-34*299792458*100/1.380649e-23
    for T in plan['comparison_T_K']:
        value=native(T);x=c2*levels/T;w=weights*np.exp(-x);p=w/w.sum();mean=p@x
        selected.append(dict(T_K=T,native_logQ_jets=value.tolist(),
            native_internal_Cv_over_kB=(value[:,1]+value[:,2]).tolist(),
            fixed_level_Q=float(w.sum()),fixed_level_DlogQ=float(mean),
            fixed_level_internal_Cv_over_kB=float(p@((x-mean)**2)),
            published_Q=float(table[T-1,1]) if T<=20000 else None))
    cp_minus_cv=table[:,5]-table[:,6]
    save('source-column-audit.json',dict(classification='Counterexample candidate',rows=20000,unique_temperatures=True,
        finite_all_columns=True,positive_Q=True,level_rows=348,unique_level_keys=True,
        source_first_row_Cp=float(table[0,5]),source_first_row_Cv=float(table[0,6]),
        min_Cp_minus_Cv=float(cp_minus_cv.min()),max_Cp_minus_Cv=float(cp_minus_cv.max()),
        tabulated_R_inferred_from_first_Cp_over_2point5=float(table[0,5]/2.5),
        warning='The source column called Cv vanishes at 1 K while Cp is 2.5 R, so it is not the total ideal-molecule Cv of article equation 45. No silent 1.5 R correction was made. High-temperature edge/derivative columns require an independent source-method audit before using them as derivative truth. Compare Q directly; retain raw columns.'))
    save('interval.json',interval_audit(rows,coefs,plan))
    w=sp.symbols('w0:3',positive=True);energies=sp.symbols('E0:3',real=True);Q=sum(w)
    lhs=Q*sum(w[i]*energies[i]**2 for i in range(3))-sum(w[i]*energies[i] for i in range(3))**2
    rhs=sum(w[i]*w[j]*(energies[i]-energies[j])**2 for i in range(3) for j in range(i+1,3))
    assert sp.expand(lhs-rhs)==0
    save('identity.json',dict(classification='Proven',passed=True,
        theorem='For Q(T)=sum g_i exp[-E_i/(k_B T)] with fixed energies, positive T-independent g_i and convergent differentiated sums, DlnQ+D^2lnQ=Var(E)/(k_B T)^2 >= 0. The pairwise numerator is sum_{i<j} w_i*w_j*(E_i-E_j)^2. It is invariant under a constant energy-zero shift. The 3-level symbolic polynomial check is a finite control of the general pairwise identity.',
        boundary='A negative internal partition Cv excludes such a positive fixed-spectrum representation. It does not prove negative total mixture heat capacity, instability of the star, or applicability when energies/occupation probabilities depend on temperature. No plasma continuum physics is certified by a finite spectral sum.'))
    comparison=np.exp(values[:,0,0])/table[999:,1]-1
    np.savez_compressed(OUT/'comparison.npz',table=table,temperatures=temperatures,native=values,replayed=replayed,
        level_term_values_cm_inverse=levels,level_weights=weights)
    save('result.json',dict(classification='Counterexample candidate',completed=True,selected=selected,
        native_replay_maximum_absolute_error=error,overlap_points=len(temperatures),
        native_to_published_Q_ratio_minus_one_min=float(comparison.min()),
        native_to_published_Q_ratio_minus_one_max=float(comparison.max()),
        stored_GR_cells_in_published_T_domain=inside.tolist(),
        stored_GR_temperature_range_K=[float(actual.min()),float(actual.max())],
        positive_fixed_spectrum_compatible_high_T_Taylor=False,physical_EOS_certified=False,
        native_EOS_replaced=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify();print('MOLECULAR PARTITION DATA',error,len(temperatures),selected[-1],flush=True)


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in json.loads((OUT/'plan.json').read_text())['runtime'].items():assert g.c.sha(path)==digest,path
    assert json.loads((OUT/'identity.json').read_text())['passed']
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS molecular source, interval model and native comparison bindings; physical EOS remains open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
