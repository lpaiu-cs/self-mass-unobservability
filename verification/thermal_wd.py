"""Request 18: source-bound thermal WD selection and linear response audit."""
from pathlib import Path
import gzip, hashlib, json, shutil, struct, subprocess, sys, urllib.request
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/thermal-wd18'
SOURCE = OUT / 'sources'
CACHE = Path('/home/lpaiu/work/thermal-wd18')
BASE = 'https://cdsarc.cds.unistra.fr/ftp/J/A+A/595/A35/'
TARGET = .197536385307


def sha(path):
    with Path(path).open('rb') as f:
        h = hashlib.sha256()
        for b in iter(lambda: f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()


def save(name, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/name).write_text(json.dumps(obj, ensure_ascii=False, indent=2)+'\n')


def fetch(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        temp = dest.with_suffix(dest.suffix+'.part')
        with urllib.request.urlopen(url, timeout=45) as r, temp.open('wb') as f:
            shutil.copyfileobj(r, f)
        temp.rename(dest)
    return dest


def acquire():
    SOURCE.mkdir(parents=True, exist_ok=True)
    plan = OUT/'preselection-plan.md'
    if not plan.exists(): shutil.copy2(ROOT/'notes/REQUEST18_THERMAL_WD_KO.md',plan)
    for source, target in [('zenodo18.tmp.json','zenodo.json'),('cds18.tmp.txt','ReadMe'),
                           ('cds18-list.tmp.txt','list.dat'),('cds18-index.tmp.html','index.html'),
                           ('cds18-rot.tmp.html','rot-index.html')]:
        if (ROOT/source).exists(): (ROOT/source).rename(SOURCE/target)
    meta = json.loads((SOURCE/'zenodo.json').read_text())
    for f in meta['files']:
        if f['size'] < 100000:
            p = fetch(f['links']['self'], SOURCE/f['key'])
            assert hashlib.md5(p.read_bytes()).hexdigest() == f['checksum'].split(':')[1]
    candidates = []
    for line in (SOURCE/'list.dat').read_text().splitlines():
        rel = line.split('|')[0].strip()
        # Filename masses have 0.001 Msun precision: include rounding margin.
        if abs(float(rel.rsplit('_',1)[1])-TARGET) <= .01*TARGET+.0005:
            fetch(BASE+rel+'.gz', SOURCE/(rel+'.gz'))
            candidates.append(rel)
    save('selection-inputs.json',dict(classification='Imported from prior work',
      catalogue_entries=266, filename_mass_prefilter=candidates,
      rule='abs(rounded_mass-target)<=0.01*target+0.0005; actual row mass tested next',
      plan_sha256=sha(plan)))
    print('Fetched',len(candidates),'candidate tracks',flush=True)


def select():
    result = []
    for rel in json.loads((OUT/'selection-inputs.json').read_text())['filename_mass_prefilter']:
        with gzip.open(SOURCE/(rel+'.gz'),'rt') as f: lines = f.readlines()
        rows = np.loadtxt(lines)
        assert rows.ndim == 2 and rows.shape[1] == (23 if rel.startswith('basic/') else 27)
        i = 6 if rel.startswith('basic/') else 9
        t, g = 10**rows[:,i], rows[:,i+3]
        mass = abs(rows[:,2]/TARGET-1) <= .01
        optical = (abs(t-15800) <= 300) & (abs(g-5.82) <= .15)
        mask = mass & optical
        distance = ((t-15800)/100)**2 + ((g-5.82)/.05)**2
        indices = np.flatnonzero(mass)
        best = int(indices[np.argmin(distance[indices])]) if len(indices) else None
        record = dict(file=rel,rows=len(rows),mass_range=[float(rows[:,2].min()),float(rows[:,2].max())],
          matching_rows=int(mask.sum()),best_index=best,columns_source='CDS ReadMe; no header row in track')
        if best is not None:
            record['best'] = dict(mass=float(rows[best,2]),Teff=float(t[best]),logg=float(g[best]),
              radius=float(10**rows[best,i+2]),age=float(rows[best,1]),model=int(rows[best,0]),
              diagnostic_chi2=float(distance[best]),raw=rows[best].tolist())
        result.append(record)
        if mask.any(): np.savetxt(OUT/(Path(rel).name+'-candidates.txt'),rows[mask],header='Columns: CDS ReadMe; original numeric rows')
    save('track-selection.json',dict(classification='Counterexample candidate',tracks=result))
    print(json.dumps([{k:v for k,v in x.items() if k!='columns'} for x in result],indent=2))


def archive_index():
    header = (SOURCE/'7z-header.bin').read_bytes()
    assert header[:6] == b'7z\xbc\xaf\x27\x1c'
    offset, size = struct.unpack('<QQ',header[12:28])
    length = 32+offset+size
    start = max(32, length-1024*1024)
    url = 'https://zenodo.org/records/2634020/files/rotation_diffusion_3.4.7z'
    req = urllib.request.Request(url, headers={'Range':f'bytes={start}-{length-1}'})
    with urllib.request.urlopen(req,timeout=45) as r:
        assert r.status == 206
        tail = r.read()
    assert len(tail) == length-start
    (SOURCE/'7z-tail.bin').write_bytes(tail)
    path = CACHE/'index-only.7z'
    with path.open('wb') as f:
        f.write(header); f.seek(start); f.write(tail)
    run = subprocess.run([str(CACHE/'p7zip/usr/lib/p7zip/7z'),'l','-slt',str(path)],capture_output=True,text=True)
    (OUT/'archive-index.txt').write_text(run.stdout+run.stderr)
    print(run.stdout[:1500],run.stderr)
    assert run.returncode == 0


def blocks():
    entries = []
    for block in (OUT/'archive-index.txt').read_text().split('\n\n'):
        d = dict(line.split(' = ',1) for line in block.splitlines() if ' = ' in line)
        if 'Block' in d: entries.append(d)
    offset = 32
    for block in range(7):
        rows = [x for x in entries if x['Block'] == str(block)]
        size = sum(int(x['Packed Size']) for x in rows if x['Packed Size'])
        print(block,offset,size,rows[0]['Path'],rows[-1]['Path'])
        for x in rows:
            if any(n in x['Path'] for n in ['profiles.index','history.data','profile379.','profile380.']): print(x['Path'])
        offset += size


def archive_fetch():
    from concurrent.futures import ThreadPoolExecutor
    url = 'https://zenodo.org/records/2634020/files/rotation_diffusion_3.4.7z'
    # ponytail: only blocks 2 and 5 contain the target-era profiles and index.
    # This sparse archive is never claimed to pass the whole-archive MD5 check.
    ranges = []
    for start, length in [(577170354,261390212),(1360912649,57047088)]:
        edges = np.linspace(start,start+length,5,dtype=np.int64)
        ranges.extend((int(a),int(b)-1) for a,b in zip(edges[:-1],edges[1:]))
    def part(pair):
        a,b = pair; path = CACHE/f'range-{a}-{b}.bin'
        if not path.exists():
            req = urllib.request.Request(url,headers={'Range':f'bytes={a}-{b}'})
            with urllib.request.urlopen(req,timeout=60) as r:
                assert r.status == 206 and r.headers['Content-Range'].startswith(f'bytes {a}-{b}/')
                with path.with_suffix('.part').open('wb') as f: shutil.copyfileobj(r,f)
            path.with_suffix('.part').rename(path)
        assert path.stat().st_size == b-a+1
        print('range complete',a,b,flush=True)
        return dict(start=a,end=b,sha256=sha(path),cache=str(path))
    with ThreadPoolExecutor(max_workers=4) as pool: records = list(pool.map(part,ranges))
    with (CACHE/'index-only.7z').open('r+b') as f:
        for r in records:
            f.seek(r['start'])
            with Path(r['cache']).open('rb') as source: shutil.copyfileobj(source,f)
    save('archive-ranges.json',dict(url=url,whole_archive_md5_checked=False,ranges=records))
    cmd = [str(CACHE/'p7zip/usr/lib/p7zip/7z'),'x','-y',str(CACHE/'index-only.7z'),
      '-o'+str(CACHE/'extracted'),'rotation_diffusion_3.4/LOGS1/profiles.index']+[
      f'rotation_diffusion_3.4/LOGS1/profile{i}.data' for i in range(314,419)]
    run = subprocess.run(cmd,capture_output=True,text=True)
    (OUT/'archive-extract.log').write_text(run.stdout+run.stderr)
    print(run.stdout[-1000:],run.stderr[-1000:])
    assert run.returncode == 0


def archive_prefix():
    path = CACHE/'rotation_diffusion_3.4.7z.part'
    assert path.stat().st_size >= 314360261
    header = (SOURCE/'7z-header.bin').read_bytes()
    offset, size = struct.unpack('<QQ',header[12:28]); length = 32+offset+size
    tail = (SOURCE/'7z-tail.bin').read_bytes()
    with path.open('r+b') as f: f.seek(length-len(tail)); f.write(tail)
    cmd = [str(CACHE/'p7zip/usr/lib/p7zip/7z'),'x','-y',str(path),'-o'+str(CACHE/'extracted'),
      'rotation_diffusion_3.4/LOGS1/history.data','rotation_diffusion_3.4/LOGS1/profile100.data']
    run = subprocess.run(cmd,capture_output=True,text=True)
    (OUT/'archive-prefix-extract.log').write_text(run.stdout+run.stderr)
    print(run.stdout[-800:],run.stderr[-800:])
    assert run.returncode == 0


def mesa(path):
    opener = (lambda:gzip.open(path,'rt')) if path.suffix == '.gz' else path.open
    with opener() as f:
        lines = [next(f) for _ in range(6)]
    header = dict(zip(lines[1].split(),lines[2].split()))
    with opener() as f: data = np.loadtxt(f,skiprows=6)
    names = lines[5].split()
    assert len(names) == data.shape[1]
    return header,{k:data[:,i] for i,k in enumerate(names)}


def inspect():
    logs = CACHE/'extracted/rotation_diffusion_3.4/LOGS1'
    index = np.loadtxt(logs/'profiles.index',skiprows=1,dtype=int)
    print('near target models:',index[abs(index[:,0]-18969)<150].tolist())
    h,d = mesa(logs/'profile100.data')
    print('header',h);print('columns',list(d))
    for model,priority,number in index[abs(index[:,0]-18969)<150]:
        path = logs/f'profile{number}.data'
        if path.exists():
            h,d = mesa(path)
            print(number,{k:h[k] for k in ['model_number','star_mass','star_age','Teff','photosphere_r']})


def profiles():
    logs = CACHE/'extracted/rotation_diffusion_3.4/LOGS1'
    shutil.copy2(logs/'profiles.index',SOURCE/'profiles.index')
    index = np.loadtxt(logs/'profiles.index',skiprows=1,dtype=int)
    with gzip.open(SOURCE/'rot/z_0.02_rotation_0.198.gz','rt') as f: track = np.loadtxt(f)
    by_model = {int(row[0]):row for row in track}
    choices = []
    for model,priority,number in index:
        if model not in by_model: continue
        r = by_model[model]; t = 10**r[9]; g = r[12]
        if abs(r[2]/TARGET-1)>.01: continue
        choices.append(dict(model=int(model),profile=int(number),mass=float(r[2]),Teff=float(t),logg=float(g),
          diagnostic_chi2=float(((t-15800)/100)**2+((g-5.82)/.05)**2),
          optical_candidate=bool(abs(t-15800)<=300 and abs(g-5.82)<=.15)))
    choices.sort(key=lambda x:x['diagnostic_chi2'])
    print('nearest saved profiles',json.dumps(choices[:10],indent=2))
    for number in [387,388]:
        src = logs/f'profile{number}.data'
        # Retain the exact published bytes, gzip-compressed locally with fixed timestamp.
        with (SOURCE/f'profile{number}.data.gz').open('wb') as f:
            with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0) as z:
                with src.open('rb') as original: shutil.copyfileobj(original,z)
    h,d = mesa(logs/'history.data')
    for model in [18950,18969,19000]:
        at=np.flatnonzero(d['model_number']==model);assert len(at)==1
        row=by_model[model]; n=at[0]
        for key,k in [('star_mass',2),('log_Teff',9),('log_R',11),('log_g',12)]:
            assert np.isclose(d[key][n],row[k],rtol=1e-12,atol=1e-12),(model,key)
    save('profile-selection.json',dict(classification='Counterexample candidate',
      saved_profile_optical_candidates=[x for x in choices if x['optical_candidate']],nearest=choices[:10],
      selected_bracket_profiles=[387,388],target_model=18969,
      history_cds_values_crosschecked=True,history_sha256=sha(logs/'history.data')))


def shell_response(r, m, omega=0., ctx=None):
    """Exact transfer through the declared constant-density Newtonian shells.

    Input r and m are geometric metres, treated as exact binary floats.
    Return static/outgoing l=0 response with incoming j0(omega*r/c).
    """
    import mpmath as mp
    from stellar_matching import C
    ctx = ctx or mp.mp
    R=ctx.mpf(float(r[-1])); wave=ctx.mpf(float(omega))/C*R
    u,v=ctx.mpf(0),ctx.mpf(1); eta=ctx.mpf(0)
    for a,b,ma,mb in zip(r[:-1],r[1:],m[:-1],m[1:]):
        a,b,ma,mb=map(lambda x:ctx.mpf(float(x)),[a,b,ma,mb])
        potential=12*(mb-ma)/(b**3-a**3) # 4*pi*|beta|*rho_geom, beta=-4.
        eta += potential*(b*b-a*a)/2
        k=ctx.sqrt(wave*wave+potential*R*R); dx=(b-a)/R
        if k == 0: u=u+dx*v
        else:
            co,si=ctx.cos(k*dx),ctx.sin(k*dx)
            u,v=co*u+si/k*v,-k*si*u+co*v
    if omega == 0: return R*(u-v)/v,eta
    f=ctx.exp(-1j*wave)*(u*ctx.cos(wave)-v*ctx.sin(wave)/wave)/(v-1j*wave*u)
    return R*f,eta


def response(nominal=False):
    import mpmath as mp
    from stellar_matching import GM_SUN,C,MSUN
    from nonzero_drive import kepler_geometry
    mp.mp.dps=60;mp.iv.dps=60
    records=[];period=kepler_geometry()[0]['period_days']*86400
    for number in [387,388]:
        path=SOURCE/f'profile{number}.data.gz';h,d=mesa(path)
        radius_unit=695700000. if nominal else 695980000.
        mass_unit=MSUN if nominal else 6.67430e-11*1.9892e30/C**2
        r=np.r_[0,d['radius'][::-1]*radius_unit]
        m=np.r_[0,d['mass'][::-1]*mass_unit]
        assert np.all(np.diff(r)>0) and np.all(np.diff(m)>0)
        assert abs(m[-1]/mass_unit-float(h['star_mass'])) < 1e-14
        np.savez_compressed(OUT/f'shells-{number}{"-nominal" if nominal else ""}.npz',r_m=r,m_geom_m=m)
        chi,eta=shell_response(r,m)
        ivchi,iveta=shell_response(r,m,ctx=mp.iv)
        assert ivchi.a <= chi <= ivchi.b and iveta.a <= eta <= iveta.b
        q0=mp.mpf(float(m[-1]))*4
        assert q0 < chi < q0/(1-eta) and eta < 1
        refinements=[]
        for stride in [2,4]:
            idx=np.unique(np.r_[np.arange(0,len(r),stride),len(r)-1])
            q,_=shell_response(r[idx],m[idx]);refinements.append(float(q/chi-1))
        dynamic=[]
        for frequency in [2*np.pi/period,2*np.pi/(327.2551609467*86400),.01,1.]:
            value,_=shell_response(r,m,frequency)
            k=mp.mpf(float(frequency))/C; kr=k*float(r[-1])
            bound=kr**2/(3*(1-eta))+abs(k)*q0/(1-eta)**2
            ik=mp.iv.mpf(float(frequency))/C;iq0=4*mp.iv.mpf(float(m[-1]));iR=mp.iv.mpf(float(r[-1]))
            ibound=(ik*iR)**2/(3*(1-iveta))+abs(ik)*iq0/(1-iveta)**2
            measured=abs(value-chi)/chi
            assert measured <= bound
            dynamic.append(dict(omega_s=frequency,response_real_m=float(value.real),response_imag_m=float(value.imag),
              relative_static_difference=float(measured),conditional_relative_bound=float(bound),
              conditional_bound_interval=str(ibound),bound_upper_rounded_outward=float(np.nextafter(float(ibound.b),np.inf))))
        shells_rho=np.diff(m)/((4*np.pi/3)*np.diff(r**3))*C*C/6.67430e-11
        published_rho=10**d['logRho'][::-1]*1000
        # Cell-averaged density diagnostic; nominal units and rounding affect it.
        density_relative=shells_rho/published_rho-1
        p=10**d['logP'][::-1]*.1;T=10**d['logT'][::-1]
        record=dict(profile=number,model=int(h['model_number']),cells=len(r)-1,
          mass_source_units=float(h['star_mass']),mass_solar=float(m[-1]/MSUN),mass_relative_target=float(m[-1]/MSUN/TARGET-1),
          Teff_K=float(h['Teff']),radius_m=float(r[-1]),radius_solar=float(h['photosphere_r']),
          central_temperature_K=float(T[0]),central_density_kg_m3=float(published_rho[0]),
          total_hydrogen_solar=float(h['star_mass_h1']),
          max_pressure_rest_energy_ratio=float(np.max(p/(published_rho*C*C))),
          max_rotation_ratio=float(np.max(abs(d['omega_div_omega_crit']))),
          density_reconstruction_relative_quantiles=np.quantile(density_relative,[0,.5,1]).tolist(),
          density_reconstruction_mass_weighted_absolute_relative=float(np.sum(abs(density_relative)*np.diff(m))/m[-1]),
          susceptibility_m=str(chi),susceptibility_interval_m=str(ivchi),interval_width_m=float(ivchi.delta),
          eta=str(eta),eta_interval=str(iveta),born_charge_m=str(q0),
          fractional_structure_correction=float(chi/q0-1),coarsening_relative_changes=refinements,
          dynamic=dynamic,input_gzip_sha256=sha(path))
        records.append(record);print(json.dumps(record,indent=2),flush=True)
    save('response-nominal.json' if nominal else 'response.json',dict(classification='Counterexample candidate',
      model='flat-space linear massless DEF scalar, beta=-4, zero background, frozen positive density shells',
      units=dict(GM_sun_SI=GM_SUN,R_source_unit_m=radius_unit,M_source_unit_geom_m=mass_unit,C_SI=C,
        adopted_G_SI=6.67430e-11,unit_reconstruction='nominal comparison' if nominal else 'source-column dimensional identities'),
      arithmetic_precision=dict(mpf_dps=mp.mp.dps,interval_dps=mp.iv.dps),
      original_MESA_source_constants_file_obtained=False,records=records))


def response_nominal(): response(True)


def calibration():
    rows=[]
    for number in [387,388]:
        h,d=mesa(SOURCE/f'profile{number}.data.gz')
        radius_unit=d['v_rot']*1000/d['omega']/d['radius']
        r=np.r_[d['radius'],0.]*float(np.median(radius_unit))
        m=np.r_[d['mass'],0.]
        mass_unit=(4*np.pi/3)*(r[:-1]**3-r[1:]**3)*(10**d['logRho']*1000)/(m[:-1]-m[1:])
        mask=(m[:-1]-m[1:])>float(h['star_mass'])*1e-7
        row=dict(profile=number,radius_unit_m_quantiles=np.quantile(radius_unit,[0,.5,1]).tolist(),
          mass_unit_kg_quantiles=np.quantile(mass_unit[mask],[0,.5,1]).tolist(),mass_cells=int(mask.sum()))
        print(row);rows.append(row)
    save('unit-diagnostics.json',dict(classification='Counterexample candidate',
      assumption='v_rot in km/s and omega in rad/s at same face; logRho is cell mean mass density',records=rows))


def check():
    import mpmath as mp, sympy as s, zlib
    from scipy.integrate import solve_ivp
    mp.mp.dps=70;mp.iv.dps=70
    k,x,u,v=s.symbols('k x u v',real=True)
    mat=s.Matrix([[s.cos(k*x),s.sin(k*x)/k],[-k*s.sin(k*x),s.cos(k*x)]])
    assert s.simplify(mat.det()-1)==0
    assert s.simplify(s.diff((mat*s.Matrix([u,v]))[0],x,2)+k*k*(mat*s.Matrix([u,v]))[0])==0
    # Independent uniform-sphere solution, vacuum, and direct ODE controls.
    rr=np.linspace(0,1.,101);mm=.001*rr**3
    value,eta=shell_response(rr,mm);z=mp.sqrt(12*mp.mpf(float(mm[-1])) )
    exact=mp.tan(z)/z-1
    assert abs(value/exact-1)<mp.mpf('1e-14')
    ode=solve_ivp(lambda t,y:[y[1],-.012*y[0]],(0,1),[0.,1.],rtol=2e-13,atol=1e-15,method='DOP853')
    assert ode.success and abs(float(value)-(ode.y[0,-1]/ode.y[1,-1]-1))<1e-13
    assert shell_response(rr,np.zeros_like(rr))[0]==0
    assert abs(shell_response(rr,np.zeros_like(rr),1e8)[0])<mp.mpf('1e-60')
    for omega in [1e4,1e8]:
        q,_=shell_response(rr,mm,omega);qm,_=shell_response(rr,mm,-omega)
        assert abs(qm-mp.conj(q))<mp.mpf('1e-60')
        assert abs(q.imag-mp.mpf(omega)/299792458*abs(q)**2)<mp.mpf('1e-60')
    selections=json.loads((OUT/'selection-inputs.json').read_text())
    assert sha(OUT/'preselection-plan.md')==selections['plan_sha256']
    source_checks=[]
    for number in [387,388]:
        path=SOURCE/f'profile{number}.data.gz'
        key=f'Path = rotation_diffusion_3.4/LOGS1/profile{number}.data\n'
        entry=(OUT/'archive-index.txt').read_text().split(key)[1].split('\n\n')[0]
        fields=dict(line.split(' = ',1) for line in entry.splitlines() if ' = ' in line)
        crc=0;size=0;digest=hashlib.sha256()
        with gzip.open(path,'rb') as f:
            for chunk in iter(lambda:f.read(1024*1024),b''):
                crc=zlib.crc32(chunk,crc);size+=len(chunk);digest.update(chunk)
        assert size==int(fields['Size']) and f'{crc:08X}'==fields['CRC']
        source_checks.append(dict(profile=number,raw_sha256=digest.hexdigest(),bytes=size,archive_crc32=fields['CRC']))
    unit=json.loads((OUT/'unit-diagnostics.json').read_text())
    for row in unit['records']:
        assert max(abs(x/695980000.-1) for x in row['radius_unit_m_quantiles'])<1e-14
        assert max(abs(x/1.9892e30-1) for x in row['mass_unit_kg_quantiles'])<2e-9
    results=json.loads((OUT/'response.json').read_text())
    for row in results['records']:
        assert row['mass_relative_target']>.002
        assert row['interval_width_m']<1e-50
        assert max(abs(x) for x in row['coarsening_relative_changes'])<1e-9
        assert row['dynamic'][0]['conditional_relative_bound']<2.1e-10
    assert not json.loads((OUT/'profile-selection.json').read_text())['saved_profile_optical_candidates']
    save('audit.json',dict(classification='Proven',uniform_sphere_analytic_and_ODE=True,
      vacuum_zero_response=True,frequency_conjugacy_and_elastic_unitarity=True,symbolic_transfer=True,
      exact_selected_source_bytes=source_checks,source_units_crosschecked=True,
      whole_archive_md5_checked=False,scope='declared shell model and source integrity, not an interval MESA evolution certificate'))
    print('PASS: 기호식, 균일 구·진공·독립 ODE·unitarity, 원자료 CRC 및 SHA, 단위와 미완료 경계')


def maintain():
    additions={
      'model-definition':'분류: Counterexample candidate. 공개 MESA 열·외피 진화 자료에서 안쪽 WD의 광학 후보 이력 행과 그 전후의 내부 구조 두 개를 확보했다. 회전속도·각속도와 셀 질량·밀도로 자료 단위를 복원했다. 두 구조의 환산 Newtonian 질량은 약 0.198040 태양질량으로 타이밍 기준보다 0.2548% 높다. 영 배경·고정 밀도·평탄 시공간의 선형 scalar 모형을 별도로 정의했으며 완전한 열 GR 별로 취급하지 않는다.',
      'observable-targets':'분류: Proven. 지정된 두 열 WD 구각 모형에서 l=0 scalar 산란 응답의 정적 값은 약 1169.86 m다. 안쪽 궤도 주파수에서 정적 값과의 상대 차이는 2.1e−10 미만이라는 조건부 해석 상계를 얻었다. 이것은 해당 독립 scalar 모형의 산란 응답이며 삼중계 힘·타이밍 잔차 또는 관측 검출값이 아니다.',
      'adiabatic-limit':'분류: Proven. 비음수·구대칭·고정 밀도에서 eta=4pi|beta|G/c² integral rho*r dr<1이면 scalar 적분 연산자가 축약 사상이어서 정적 해가 유일하다. Q0=|beta|GM/c²에 대해 Q0<=chi<=Q0/(1−eta), |f(k)−chi|/chi<=(kR)²/[3(1−eta)]+|k|Q0/(1−eta)²를 얻었다. 두 공개 열 구조에서 정의한 구각 모형은 eta<0.000166이다. 실제 시간 의존 항성 구조를 이 두 끝점으로 포괄한 것은 아니다.',
      'nonadiabatic-regime':'분류: Proven. 고정 밀도의 약한 중력 열 WD scalar 모형에서 궤도 주파수 응답은 정적 값의 2.1e−10 이내이므로 큰 비단열 응답을 공급하지 않는다. 상반평면의 독립 scalar 성장 해도 축약 조건으로 배제된다. 하반평면의 모든 pole이나 유체·회전·중성자별 모드를 배제하는 정리는 아니다. 분류: Conjectural. 비영 배경의 유체·metric 결합과 실제 상호 구동은 추가 matching이 필요하다.',
      'failure-ledger-dynamic-chi':'분류: Proven. 공개 열 진화 이력에는 광학 후보가 있지만 해당 실행의 저장된 내부 구조 중 사전 질량·광학 기준을 함께 만족하는 것은 0개였다. 최적 이력 행의 전후 구조는 온도 16409 K와 14577 K로 각각 기준을 벗어난다. 단위 보정 뒤 질량도 목표보다 0.2548% 높다. 임의 보간이나 질량 재규격화를 실제 EOS matching으로 승격하지 않았다. 조건부 고정 구각 모형의 산술 구간·주파수 상계는 확보했으나 실제 열 별의 구간 인증과 물리 timing·전체 비선형 추론은 미완료다.'}
    old=json.loads((ROOT/'outputs/nonzero-drive17/manifest.json').read_text())['sha256']
    rels=['docs/'+x+'.md' for x in additions]+['paper/revision-manifest.json']
    dest=OUT/'request17-notes';dest.mkdir(exist_ok=False);bindings={}
    for rel in rels:
        assert sha(ROOT/rel)==old[rel]
        snap=dest/Path(rel).name;shutil.copy2(ROOT/rel,snap)
        bindings[rel]=dict(snapshot=snap.relative_to(ROOT).as_posix(),sha256=old[rel],historical_manifest=rel.startswith('paper/'))
    save('historical-note-bindings.json',bindings)
    revision=json.loads((ROOT/'paper/revision-manifest.json').read_text())
    for stem,body in additions.items():
        rel='docs/'+stem+'.md'
        with (ROOT/rel).open('ab') as f:f.write(('\n\n## Request 18 열 백색왜성 구조와 응답 경계\n\n'+body+'\n\n세부 근거: [한글 도출·검증 보고서](../notes/REQUEST18_THERMAL_WD_KO.md).\n').encode())
        revision['sha256'][rel]=sha(ROOT/rel)
    revision['request18_supporting_note_update']=dict(evidence_manifest='outputs/thermal-wd18/manifest.json',
      historical_notes='outputs/thermal-wd18/historical-note-bindings.json',
      status='공개 열 항성 구조 확보와 고정 구각 scalar 응답의 조건부 정리; 정확한 질량·광학 구조와 전체 timing·추론 미완료',
      artifact_status='원고 PDF와 ZIP은 Request 12 동결본이다.')
    (ROOT/'paper/revision-manifest.json').write_text(json.dumps(revision,ensure_ascii=False,indent=2)+'\n')


def seal():
    import scipy,mpmath,sympy,nonzero_drive,stellar_matching,remaining_audit
    modules=[np,scipy,mpmath,sympy,nonzero_drive,stellar_matching,remaining_audit]
    save('provenance.json',dict(before_task_checkpoint='5f99017',interpreter=sys.executable,
      dependency_path='/home/lpaiu/work/nutimo_pilot/request13_deps',
      versions={m.__name__:m.__version__ for m in [np,scipy,mpmath,sympy]},
      module_file_sha256={str(m.__file__):sha(m.__file__) for m in modules},
      input_sha256={'outputs/validated-variational/ivp.hex':sha(ROOT/'outputs/validated-variational/ivp.hex')},
      producer='verification/thermal_wd.py select; profiles; calibration; response_nominal; response; check'))
    save('gates.json',dict(classification='Proven',theorem_progress=True,published_thermal_profiles_acquired=True,
      source_dimensions_crosschecked=True,declared_shell_static_arithmetic_certified=True,
      declared_shell_frequency_bound_certified=True,exact_J0337_mass_and_optical_profile_matched=False,
      thermal_evolution_rerun=False,full_stellar_interval_certificate=False,
      genuine_orbital_timescale_state_established=False,full_physical_dynamic_force_and_readout=False,
      full_span_variational_certificate=False,full_28_parameter_initialization_certificate=False,
      complete_nonlinear_observational_inference=False))
    paths=[p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    paths += [ROOT/'verification/thermal_wd.py',ROOT/'notes/REQUEST18_THERMAL_WD_KO.md']
    paths += [ROOT/x for x in json.loads((OUT/'historical-note-bindings.json').read_text())]
    save('manifest.json',dict(classification='Proven',scope='Request 18 열 구조 자료 및 조건부 정리',
      sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(paths)}))


def verify():
    histories={
      'outputs/validated-variational/manifest.json':'outputs/remaining-levers15/historical-note-bindings.json',
      'outputs/remaining-levers15/manifest.json':'outputs/nbody-readout16/historical-note-bindings.json',
      'outputs/nbody-readout16/manifest.json':'outputs/nonzero-drive17/historical-note-bindings.json',
      'outputs/nonzero-drive17/manifest.json':'outputs/thermal-wd18/historical-note-bindings.json'}
    count=0
    for label in ['outputs/research-remediation/manifest.json',*histories,'paper/revision-manifest.json','outputs/thermal-wd18/manifest.json']:
        old=json.loads((ROOT/histories[label]).read_text()) if label in histories else {}
        for name,expected in json.loads((ROOT/label).read_text())['sha256'].items():
            path=ROOT/name
            if name in old:
                bind=old[name];path=ROOT/bind['snapshot'];assert expected==bind['sha256']
                if bind.get('historical_manifest'):
                    before=json.loads(path.read_text());after=json.loads((ROOT/name).read_text())
                    for k,value in before.items():
                        if k!='sha256':assert after[k]==value,k
                    for k,value in before['sha256'].items():
                        if k not in old:assert after['sha256'][k]==value,k
                else:assert (ROOT/name).read_bytes().startswith(path.read_bytes())
            assert sha(path)==expected,(label,name)
            count+=1
    for name,expected in json.loads((OUT/'provenance.json').read_text())['input_sha256'].items():assert sha(ROOT/name)==expected
    gates=json.loads((OUT/'gates.json').read_text());assert gates['theorem_progress'] and not gates['complete_nonlinear_observational_inference']
    print('PASS:',count,'현재·역사적 SHA, 원고 동결 및 미완료 판정 보존')


if __name__ == '__main__': globals()[sys.argv[1]]()
