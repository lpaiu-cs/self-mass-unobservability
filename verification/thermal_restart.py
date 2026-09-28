"""Request 19: reproduce the published thermal MESA evolution before calibration."""
from pathlib import Path
import gzip,hashlib,json,os,shutil,struct,subprocess,sys,time,urllib.request,zipfile,zlib

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/thermal-restart19'
CACHE=Path('/home/lpaiu/work/thermal-restart19')
PREV=Path('/home/lpaiu/work/thermal-wd18')
SEVEN=str(PREV/'p7zip/usr/lib/p7zip/7z')


def save(name,data):
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/name).write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')


def fetch(url,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    if not path.exists():
        with urllib.request.urlopen(url,timeout=60) as r,path.with_suffix(path.suffix+'.part').open('wb') as f:
            shutil.copyfileobj(r,f)
        path.with_suffix(path.suffix+'.part').rename(path)
    return path


def acquire():
    OUT.mkdir(parents=True,exist_ok=True);CACHE.mkdir(parents=True,exist_ok=True)
    fetch('https://zenodo.org/api/records/2630796',OUT/'mesa-release.json')
    print([(f['key'],f['size'],f['checksum']) for f in json.loads((OUT/'mesa-release.json').read_text())['files']])
    url='https://zenodo.org/records/2634020/files/rotation_diffusion_3.4.7z'
    start,end=1417959737,1420084785
    request=urllib.request.Request(url,headers={'Range':f'bytes={start}-{end}'})
    with urllib.request.urlopen(request,timeout=60) as r:
        assert r.status==206
        data=r.read()
    assert len(data)==end-start+1
    # Copy the frozen source's sparse index; fill the published executable block here.
    archive=CACHE/'executable-block.7z'
    with archive.open('wb') as f:
        f.write((ROOT/'outputs/thermal-wd18/sources/7z-header.bin').read_bytes())
        tail=(ROOT/'outputs/thermal-wd18/sources/7z-tail.bin').read_bytes()
        f.seek(1420094073-len(tail));f.write(tail);f.seek(start);f.write(data)
    run=subprocess.run([SEVEN,'x','-y',str(archive),'-o'+str(CACHE/'published'),'rotation_diffusion_3.4/binary'],capture_output=True,text=True)
    (OUT/'extract-executable.log').write_text(run.stdout+run.stderr);assert run.returncode==0
    print(run.stdout[-300:])


def mesa_index():
    url='https://downloads.sourceforge.net/project/mesa/releases/mesa-r7624.zip'
    length=1065123758;start=length-4*1024*1024
    req=urllib.request.Request(url,headers={'Range':f'bytes={start}-{length-1}'})
    with urllib.request.urlopen(req,timeout=60) as response:
        assert response.status==206
        tail=response.read()
    assert len(tail)==length-start
    path=CACHE/'mesa-index.zip'
    with path.open('wb') as f:f.seek(start);f.write(tail)
    with zipfile.ZipFile(path) as z:
        infos=z.infolist()
        groups={}
        for f in infos:
            key='/'.join(f.filename.split('/')[:2]);groups[key]=groups.get(key,0)+f.compress_size
        save('mesa-zip-index.json',[dict(name=f.filename,offset=f.header_offset,size=f.file_size,compressed=f.compress_size,crc=f.CRC) for f in infos])
        print(groups)
    (CACHE/'mesa-zip-tail.bin').write_bytes(tail)


def mesa_sources():
    wanted=['const/public/const_def.f90','binary/private/binary_photos.f90','binary/private/run_binary_support.f90',
      'binary/job/run_binary.f','star/private/photo_in.f90','star/private/photo_out.f90',
      'star/defaults/controls.defaults','binary/defaults/binary_controls.defaults','binary/defaults/binary_job.defaults',
      'star/defaults/star_job.defaults','star/public/star_def.f90','utils/install_mesa','install','build_data_and_export']
    path=CACHE/'mesa-index.zip'
    infos=json.loads((OUT/'mesa-zip-index.json').read_text())
    for f in infos:
        if any(f['name'].endswith('/'+x) for x in wanted) or ('binary/private/' in f['name'] and f['name'].endswith(('.f','.f90'))):
            if (OUT/'sources'/f['name']).exists():continue
            start=f['offset'];end=min(1065123758-1,start+f['compressed']+65536)
            req=urllib.request.Request('https://downloads.sourceforge.net/project/mesa/releases/mesa-r7624.zip',headers={'Range':f'bytes={start}-{end}'})
            with urllib.request.urlopen(req,timeout=60) as response:
                assert response.status==206;data=response.read()
            assert len(data)==end-start+1
            with path.open('r+b') as out:out.seek(start);out.write(data)
            with zipfile.ZipFile(path) as z:
                data=z.read(f['name'])
                dest=OUT/'sources'/f['name'];dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)
            print('source',f['name'],flush=True)


def download(meta_path):
    from concurrent.futures import ThreadPoolExecutor
    f=next(x for x in json.loads(meta_path.read_text())['files'] if x['size']>1000000)
    length=f['size'];url=f['links']['self'].replace('/api/records/','/records/').removesuffix('/content');dest=CACHE/f['key']
    if f['key']=='mesa-r7624.zip':url='https://downloads.sourceforge.net/project/mesa/releases/mesa-r7624.zip'
    def chunk(n):
        start=length*n//8;end=length*(n+1)//8-1;path=CACHE/(f['key']+f'.chunk{n}')
        if not path.exists():
            req=urllib.request.Request(url,headers={'Range':f'bytes={start}-{end}'})
            with urllib.request.urlopen(req,timeout=90) as r:
                assert r.status==206
                with path.with_suffix(path.suffix+'.part').open('wb') as out:shutil.copyfileobj(r,out)
            path.with_suffix(path.suffix+'.part').rename(path)
        assert path.stat().st_size==end-start+1
        print('chunk',f['key'],n,'complete',flush=True)
        return path
    if not dest.exists():
        with ThreadPoolExecutor(max_workers=4) as pool:parts=list(pool.map(chunk,range(8)))
        with dest.open('wb') as out:
            for p in parts:
                with p.open('rb') as inp:shutil.copyfileobj(inp,out)
    digest=hashlib.md5()
    with dest.open('rb') as inp:
        for block in iter(lambda:inp.read(1024*1024),b''):digest.update(block)
    assert digest.hexdigest()==f['checksum'].split(':')[1]
    print('MD5 verified',dest,flush=True)


def sdk_download():download(OUT/'sdk-record.json')
def mesa_download():download(OUT/'mesa-release.json')


def restore_inputs():
    published=CACHE/'published/rotation_diffusion_3.4'
    for item in json.loads((OUT/'source-bindings.json').read_text())['archive_members']:
        rel=Path(item['archive_member']).relative_to('rotation_diffusion_3.4')
        dest=(published/rel).resolve();assert dest.is_relative_to(published.resolve())
        src=ROOT/item['snapshot']
        raw=gzip.decompress(src.read_bytes()) if src.suffix=='.gz' else src.read_bytes()
        assert hashlib.sha256(raw).hexdigest()==item['raw_sha256'] and f'{zlib.crc32(raw):08X}'==item['crc32']
        if dest.exists():assert dest.read_bytes()==raw
        else:dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(raw)
    (published/'binary').chmod(0o755)
    print('PASS: author input snapshots restored and checked')


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()


def prepare():
    photo=sys.argv[2];stop=int(sys.argv[3]);name=sys.argv[4]
    assert photo in ('18000','19000') and name.replace('_','').isalnum()
    assert int(photo)<=stop<=19500
    published=CACHE/'published/rotation_diffusion_3.4';run=CACHE/name
    run.mkdir(exist_ok=False)
    for p in published.iterdir():
        if p.is_file() and (p.name.startswith('inlist') or p.suffix=='.list' or p.name=='binary'):
            shutil.copy2(p,run/p.name)
    (run/'photos1').mkdir();shutil.copy2(published/'photos1'/photo,run/'photos1'/photo)
    (run/'.restart').write_text(photo+'\n')
    p=run/'inlist1';original=p.read_text()
    assert original.count('max_model_number = 100000')==1
    p.write_text(original.replace('max_model_number = 100000',f'max_model_number = {stop}').replace('profile_interval = 50','profile_interval = 1'))
    save(name+'-inputs.json',dict(classification='Counterexample candidate',photo=photo,stop=stop,
      run_directory=str(run),changes=['max_model_number','profile_interval'],
      sha256={str(p.relative_to(run)):sha(p) for p in sorted(run.rglob('*')) if p.is_file()}))
    plan=OUT/'preregistered-plan.md'
    if not plan.exists():shutil.copy2(ROOT/'notes/REQUEST19_THERMAL_RESTART_KO.md',plan)
    print(run)


def bind_sources():
    published=CACHE/'published/rotation_diffusion_3.4'
    entries={}
    for block in (ROOT/'outputs/thermal-wd18/archive-index.txt').read_text().split('\n\n'):
        fields=dict(line.split(' = ',1) for line in block.splitlines() if ' = ' in line)
        if 'Path' in fields:entries[fields['Path']]=fields
    bound=[]
    for p in sorted(published.rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(published).as_posix();entry=entries['rotation_diffusion_3.4/'+rel]
        raw=p.read_bytes();crc=f'{zlib.crc32(raw):08X}'
        assert len(raw)==int(entry['Size']) and crc==entry['CRC'],rel
        dest=OUT/'published-inputs'/rel;dest.parent.mkdir(parents=True,exist_ok=True)
        if len(raw)>100000:
            dest=dest.with_suffix(dest.suffix+'.gz')
            with dest.open('wb') as f:
                with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0) as g:g.write(raw)
        else:dest.write_bytes(raw)
        bound.append(dict(archive_member='rotation_diffusion_3.4/'+rel,bytes=len(raw),crc32=crc,
          raw_sha256=hashlib.sha256(raw).hexdigest(),snapshot=str(dest.relative_to(ROOT))))
    source=OUT/'sources/mesa-r7624/const/public/const_def.f90'
    for exact in ['standard_cgrav = 6.67428d-8','msol = 1.9892d33','rsol = 6.9598d10','secyer = 3.1558149984d7']:
        assert exact in source.read_text()
    G=6.67428e-11;M=1.9892e30;GM_sun=1.3271244e20;target=.197536385307
    source_mass=.19796106308724848
    save('source-bindings.json',dict(classification='Proven',archive_members=bound,
      whole_stellar_archive_md5_checked=False,constant_source_sha256=sha(source),
      constants=dict(G_SI=G,Msun_kg=M,Rsun_m=6.9598e8,year_s=3.1558149984e7),
      unit_comparison=dict(source_mass=source_mass,source_GM_over_adopted_GM_sun=source_mass*G*M/GM_sun,
        target_in_source_units=target*GM_sun/(G*M),relative_difference=source_mass*G*M/GM_sun/target-1),
      mass_interpretation='Newtonian model mass mapped by source G; no claim of relativistic ADM equivalence'))
    print('PASS: published archive member CRCs and original MESA constants')


def reference():
    import numpy as np
    from thermal_wd import mesa
    logs=PREV/'extracted/rotation_diffusion_3.4/LOGS1'
    assert sha(logs/'history.data')==json.loads((ROOT/'outputs/thermal-wd18/profile-selection.json').read_text())['history_sha256']
    _,d=mesa(logs/'history.data')
    keys=['model_number','star_age','star_mass','log_dt','num_zones','log_Teff','log_R',
      'log_g','log_center_T','log_center_Rho','log_center_P','total_mass_h1','period_days']
    rows={}
    for model in [18000,18950,18969,19000,19001]:
        at=np.flatnonzero(d['model_number']==model);assert len(at)==1
        rows[str(model)]={k:float(d[k][at[0]]) for k in keys}
    h,_=mesa(ROOT/'outputs/thermal-wd18/sources/profile388.data.gz')
    save('published-reference.json',dict(classification='Imported from prior work',history=rows,profile388_header=h))
    for model,row in rows.items():print(model,{k:row[k] for k in ['star_age','star_mass','log_dt','log_Teff','log_R']})


def run():
    name=sys.argv[2]
    assert name.replace('_','').isalnum()
    record=json.loads((OUT/(name+'-inputs.json')).read_text());folder=Path(record['run_directory'])
    assert folder.parent==CACHE and not (folder/'execution.log').exists()
    for rel,digest in record['sha256'].items():assert sha(folder/rel)==digest,rel
    p=folder/'inlist1'
    p.write_text(p.read_text().replace('&star_job','&star_job\n      profile_starting_model = .true.',1))
    record['changes'].append('profile_starting_model');record['sha256']['inlist1']=sha(p)
    save(name+'-inputs.json',record)
    threads=sys.argv[3] if len(sys.argv)>3 else '1';assert threads in ['1','4']
    env=os.environ.copy();env.update(MESA_DIR=str(CACHE/'mesa-r7624'),OMP_NUM_THREADS=threads,
      LD_LIBRARY_PATH=str(CACHE/'mesasdk/lib')+':'+str(CACHE/'mesasdk/lib64'))
    start=time.monotonic()
    with (folder/'execution.log').open('w') as log:
        result=subprocess.run(['./binary'],cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT)
    shutil.copy2(folder/'execution.log',OUT/(name+'-execution.log'))
    save(name+'-execution.json',dict(classification='Proven',returncode=result.returncode,
      elapsed_s=time.monotonic()-start,environment={k:env[k] for k in ['MESA_DIR','OMP_NUM_THREADS','LD_LIBRARY_PATH']},
      warning='Fortran STOP may exit zero on failure; output comparisons determine success'))
    print(name,'exit',result.returncode,'seconds',time.monotonic()-start)


def compare():
    import numpy as np
    from thermal_wd import mesa
    name=sys.argv[2];folder=CACHE/name
    original_h,original_d=mesa(ROOT/'outputs/thermal-wd18/sources/profile388.data.gz')
    restored=[]
    for path in sorted((folder/'LOGS1').glob('profile*.data')):
        try:h,d=mesa(path)
        except ValueError as e:
            restored.append(dict(file=path.name,passed=False,parse_failure=str(e)))
            continue
        if int(h['model_number'])!=19000:continue
        errors={k:abs(float(h[k])/float(original_h[k])-1) for k in ['star_mass','photosphere_r','Teff']}
        for key in ['logT','logRho']:
            errors['center_'+key]=abs(10**(d[key][-1]-original_d[key][-1])-1)
        same_mesh=len(d['logT'])==len(original_d['logT'])
        restored.append(dict(file=path.name,model=19000,relative_errors=errors,
          same_mesh=same_mesh,passed=all(v<=1e-10 for v in errors.values()),
          max_shell_logT_relative_error=float(np.max(abs(10**(d['logT']-original_d['logT'])-1))) if same_mesh else None))
    evolved=[]
    history=folder/'LOGS1/history.data'
    if history.exists():
        with history.open() as f:headers=[next(f) for _ in range(6)]
        data=np.loadtxt(history,skiprows=6,ndmin=2)
        d={k:data[:,i] for i,k in enumerate(headers[5].split())}
        expected=json.loads((OUT/'published-reference.json').read_text())['history']
        for model in [19001]:
            at=np.flatnonzero(d['model_number']==model)
            if not len(at):continue
            i=at[-1];ref=expected[str(model)]
            errors={k:abs(float(d[k][i])/ref[k]-1) for k in ['star_mass']}
            for k in ['log_R','log_Teff','log_center_T','log_center_Rho']:
                errors[k]=abs(float(10**(d[k][i]-ref[k])-1))
            dt=abs(float(10**(d['log_dt'][i]-ref['log_dt'])-1))
            evolved.append(dict(model=model,relative_errors=errors,dt_relative_error=dt,
              passed=all(v<=1e-5 for v in errors.values()) and dt<=.01))
    save(name+'-comparison.json',dict(classification='Proven',restoration=restored,one_step=evolved,
      restoration_passed=bool(restored) and all(x['passed'] for x in restored),
      one_step_passed=bool(evolved) and all(x['passed'] for x in evolved)))
    print(json.dumps(dict(restoration=restored,one_step=evolved),indent=2))


def environment():
    archives=[]
    for filename in ['mesa-release.json','sdk-record.json']:
        meta=json.loads((OUT/filename).read_text())
        item=next(x for x in meta['files'] if x['size']>1000000);p=CACHE/item['key']
        h=hashlib.md5()
        with p.open('rb') as f:
            for block in iter(lambda:f.read(1048576),b''):h.update(block)
        assert h.hexdigest()==item['checksum'].split(':')[1] and p.stat().st_size==item['size']
        archives.append(dict(name=p.name,bytes=p.stat().st_size,md5=h.hexdigest(),sha256=sha(p)))
    with zipfile.ZipFile(CACHE/'mesa-r7624.zip') as z:
        for p in (OUT/'sources').rglob('*'):
            if p.is_file():assert p.read_bytes()==z.read(p.relative_to(OUT/'sources').as_posix()),str(p)
    libraries=[]
    for line in (OUT/'published-executable-ldd.txt').read_text().splitlines():
        assert 'not found' not in line
        words=line.split()
        paths=[Path(x) for x in words if x.startswith('/')]
        libraries.extend(dict(path=str(p),sha256=sha(p)) for p in paths)
    mesa=CACHE/'mesa-r7624'
    tables={p.relative_to(mesa).as_posix():sha(p) for p in sorted((mesa/'data').rglob('*'))
      if p.is_file() and 'cache' not in p.parts}
    save('environment.json',dict(classification='Proven',archives=archives,dynamic_libraries=libraries,
      data_sha256=tables,distributed_version_number=(mesa/'data/version_number').read_text().strip(),
      release_identity='mesa-r7624 ZIP verified against Zenodo MD5; internal version_number retained as distributed',
      original_source_fragments_match_complete_archive=True))
    print('PASS: archives, source fragments, dynamic libraries,',len(tables),'data files')


def network_diagnosis():
    from thermal_wd import mesa
    p=CACHE/'published/rotation_diffusion_3.4/photos1/19000'
    records=[]
    with p.open('rb') as f:
        while marker:=f.read(4):
            n,=struct.unpack('<I',marker);payload=f.read(n)
            assert f.read(4)==marker and len(payload)==n
            records.append(payload)
    scalar=records[2]
    species,reactions,model=struct.unpack_from('<iii',scalar,60)
    assert species==22 and model==19000
    h,d=mesa(ROOT/'outputs/thermal-wd18/sources/profile388.data.gz')
    isotopes=[k for k in d if k.rstrip('0123456789') in ['h','he','c','n','o','f','ne','mg','ca'] and k[-1].isdigit()]
    save('network-diagnosis.json',dict(classification='Proven',photo_version=struct.unpack('<i',records[0])[0],
      photo_species=species,photo_reactions=reactions,photo_model=model,profile_isotopes=isotopes,
      calcium40_min=float(d['ca40'].min()),calcium40_max=float(d['ca40'].max()),
      first_run_failed=True,error='xa_old / xa_older dimensions 22 versus 21 after change_net',
      missing_original_network_file=True))
    print(species,reactions,model,isotopes)


def calcium_network():
    name=sys.argv[2];record=json.loads((OUT/(name+'-inputs.json')).read_text());folder=Path(record['run_directory'])
    assert folder.parent==CACHE and not (folder/'execution.log').exists()
    p=folder/'cno_extras.net'
    p.write_text("include 'basic.net'\ninclude 'add_cno_extras'\nadd_isos(ca40)\n")
    record['changes'].append('cno_extras.net: explicit candidate reconstruction adds inert ca40')
    record['sha256'][p.name]=sha(p);save(name+'-inputs.json',record)


def resume_controls():
    name=sys.argv[2];record=json.loads((OUT/(name+'-inputs.json')).read_text());folder=Path(record['run_directory'])
    assert folder.parent==CACHE and not (folder/'execution.log').exists()
    p=folder/'inlist1';s=p.read_text();assert s.count('do_element_diffusion = .false.')==1
    s=s.replace('do_element_diffusion = .false.','do_element_diffusion = .true.')
    assert 'which_atm_option' not in s
    p.write_text(s.replace('&controls',"&controls\n      which_atm_option = 'WD_tau_25_tables'",1))
    record['changes']+=['resume diffusion already enabled by extras_finish_step',
      'resume WD_tau_25_tables retained after earlier Teff<10000 crossing; source-derived reconstruction']
    record['sha256'][p.name]=sha(p);save(name+'-inputs.json',record)


def target_profile():
    name=sys.argv[2];record=json.loads((OUT/(name+'-inputs.json')).read_text());folder=Path(record['run_directory'])
    assert folder.parent==CACHE and not (folder/'execution.log').exists()
    p=folder/'inlist1';s=p.read_text().replace('profile_interval = 1','profile_interval = 50')
    p.write_text(s.replace('&star_job','&star_job\n      profile_model_number = 18969',1))
    record['changes']+=['retain original profile cadence 50 and save specified target model 18969']
    record['sha256'][p.name]=sha(p);save(name+'-inputs.json',record)


def trajectory():
    import numpy as np
    from thermal_wd import mesa
    _,ref=mesa(PREV/'extracted/rotation_diffusion_3.4/LOGS1/history.data')
    path=CACHE/'replay18000/LOGS1/history.data';lines=path.read_text().splitlines()
    names=lines[5].split();rows=lines[6:]
    if rows and len(rows[-1].split())!=len(names):rows=rows[:-1] # final buffered row may still be in flight
    data=np.loadtxt(rows,ndmin=2);d={k:data[:,i] for i,k in enumerate(names)}
    assert np.array_equal(d['model_number'],np.arange(18001,int(d['model_number'][-1])+1))
    lookup={int(m):i for i,m in enumerate(ref['model_number'])};ri=np.array([lookup[int(m)] for m in d['model_number']])
    errors={}
    for k in ['star_mass','log_R','log_Teff','log_center_T','log_center_Rho','log_dt']:
        v=abs(10**(d[k]-ref[k][ri])-1) if k.startswith('log_') else abs(d[k]/ref[k][ri]-1)
        errors[k]=dict(max_relative_error=float(v.max()),worst_model=int(d['model_number'][np.argmax(v)]))
    complete=int(d['model_number'][-1])==19000
    result=dict(classification='Proven',last_model=int(d['model_number'][-1]),rows=len(data),complete=complete,
      errors=errors,passed=complete and all(v['max_relative_error']<=(.01 if k=='log_dt' else 1e-5) for k,v in errors.items()))
    save('replay-comparison.json',result);print(json.dumps(result,indent=2))


def restoration_audit():
    import numpy as np
    from thermal_wd import mesa
    h,a=mesa(ROOT/'outputs/thermal-wd18/sources/profile388.data.gz')
    restored,b=mesa(CACHE/'state19000/LOGS1/profile1.data')
    keys=['logT','logRho','logP','logR','mass','radius','omega']
    keys+=json.loads((OUT/'network-diagnosis.json').read_text())['profile_isotopes']
    errors={k:float(np.max(abs(a[k]-b[k]))) for k in keys}
    assert all(np.array_equal(a[k],b[k]) for k in keys)
    assert json.loads((OUT/'state19000-comparison.json').read_text())['one_step_passed']
    replay=json.loads((OUT/'replay-comparison.json').read_text())
    if replay['complete']:
        index=np.loadtxt(CACHE/'replay18000/LOGS1/profiles.index',skiprows=1,dtype=int)
        at=index[index[:,0]==19000];assert len(at)==1
        _,final=mesa(CACHE/'replay18000/LOGS1'/f'profile{at[0,2]}.data')
        assert all(np.array_equal(a[k],final[k]) for k in keys)
    save('restoration-audit.json',dict(classification='Proven',cells=len(a['logT']),
      column_absolute_errors=errors,all_selected_structure_and_composition_columns_identical=True,
      replay_final_structure_and_composition_columns_identical=replay['complete'],
      original_hidden_network_recovered_as_bytes=False,
      conclusion='explicit reconstructed network and runtime controls reproduce saved state and next step'))
    print('PASS:',len(a['logT']),'cells,',len(keys),'structural/composition columns and one-step positive control')


def target_response():
    import numpy as np
    from thermal_wd import mesa
    assert json.loads((OUT/'replay-comparison.json').read_text())['passed']
    logs=CACHE/'replay18000/LOGS1'
    index=np.loadtxt(logs/'profiles.index',skiprows=1,dtype=int)
    at=index[index[:,0]==18969];assert len(at)==1
    path=logs/f'profile{at[0,2]}.data';h,d=mesa(path)
    ref=json.loads((OUT/'published-reference.json').read_text())['history']['18969']
    assert float(h['star_mass'])==ref['star_mass']
    assert abs(float(h['Teff'])/10**ref['log_Teff']-1)<1e-14
    assert abs(float(h['photosphere_r'])/10**ref['log_R']-1)<1e-14
    dest=OUT/'target-model18969.data.gz'
    with dest.open('wb') as f:
        with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0) as g:
            with path.open('rb') as p:shutil.copyfileobj(p,g)
    response_profile(path,ref['log_g'],'target',True)


def response_profile(path,logg,stem,optical):
    import numpy as np,mpmath as mp,sympy as sp
    from thermal_wd import mesa,shell_response
    from stellar_matching import C,GM_SUN
    from nonzero_drive import kepler_geometry
    h,d=mesa(path)
    const=json.loads((OUT/'source-bindings.json').read_text())['constants']
    unit=const['G_SI']*const['Msun_kg']/C**2
    r=np.r_[0,d['radius'][::-1]*const['Rsun_m']];m=np.r_[0,d['mass'][::-1]*unit]
    assert np.all(np.diff(r)>0) and np.all(np.diff(m)>0)
    np.savez_compressed(OUT/(stem+'-shells.npz'),r_m=r,m_geom_m=m)
    mp.mp.dps=65;mp.iv.dps=65
    chi,eta=shell_response(r,m);ivchi,iveta=shell_response(r,m,ctx=mp.iv)
    assert ivchi.a<=chi<=ivchi.b and iveta.a<=eta<=iveta.b and eta<1
    omega=2*np.pi/(kepler_geometry()[0]['period_days']*86400)
    dynamic,_=shell_response(r,m,omega)
    ik=mp.iv.mpf(float(omega))/C;iq0=4*mp.iv.mpf(float(m[-1]))
    bound=(ik*float(r[-1]))**2/(3*(1-iveta))+abs(ik)*iq0/(1-iveta)**2
    assert abs(dynamic/chi-1)<=bound.a
    # Exact algebra for the homogeneous shell propagator and unit conversion.
    z,t=sp.symbols('z t',real=True)
    transfer=sp.Matrix([[sp.cos(z*t),sp.sin(z*t)/z],[-z*sp.sin(z*t),sp.cos(z*t)]])
    assert sp.simplify(transfer.det()-1)==0
    assert sp.simplify(sp.diff(transfer[0,1],t,2)+z*z*transfer[0,1])==0
    result=dict(classification='Proven',conditional_model='frozen positive density, flat-space massless DEF beta=-4, zero background',
      model=int(h['model_number']),raw_profile_sha256=sha(path),cells=len(r)-1,Teff_K=float(h['Teff']),logg=logg,
      radius_m=float(r[-1]),mass_source_units=float(h['star_mass']),
      mass_in_GM_sun=float(m[-1]*C**2/GM_SUN),target_mass_relative_difference=float(m[-1]*C**2/GM_SUN/.197536385307-1),
      susceptibility_m=str(chi),susceptibility_interval_m=str(ivchi),static_interval_width_m=float(ivchi.delta),
      eta=str(eta),omega_s=float(omega),relative_static_difference=float(abs(dynamic/chi-1)),
      relative_bound_interval=str(bound),relative_bound_upper=float(np.nextafter(float(bound.b),np.inf)),
      optical_target_recovered_by_actual_evolution=optical,exact_mass_matched=False,
      source_GM_target_matched=bool(abs(m[-1]*C**2/GM_SUN/.197536385307-1)<=1e-8),
      full_stellar_or_gr_error_certificate=False,symbolic_propagator_checked=True)
    save(stem+'-response.json',result);print(json.dumps(result,indent=2))


def mass_response():
    result=json.loads((OUT/'mass-selection.json').read_text())
    assert result['run_completed'] and result['normal_termination']
    response_profile(Path(result['best_profile']),result['best']['logg'],'mass',result['best']['optical_candidate'])


def mass_trial():
    assert json.loads((OUT/'replay-comparison.json').read_text())['passed']
    name=sys.argv[2];record=json.loads((OUT/(name+'-inputs.json')).read_text());folder=Path(record['run_directory'])
    assert folder.parent==CACHE and not (folder/'execution.log').exists()
    target=json.loads((OUT/'source-bindings.json').read_text())['unit_comparison']['target_in_source_units']
    p=folder/'inlist1';s=p.read_text()
    p.write_text(s.replace('&star_job',f'&star_job\n      relax_mass = .true.\n      new_mass = {target:.17e}\n      lg_max_abs_mdot = -6',1))
    record['changes']+=['relax_mass via removal of envelope; rate capped at 1e-6 source Msun/year']
    record['mass_trial_plan_sha256']=sha(OUT/'mass-trial-plan.md')
    record['sha256'][p.name]=sha(p);save(name+'-inputs.json',record)


def mass_selection():
    import numpy as np
    name='mass18000';folder=CACHE/name
    lines=(folder/'LOGS1/history.data').read_text().splitlines();names=lines[5].split();rows=lines[6:]
    if len(rows[-1].split())!=len(names):rows=rows[:-1]
    data=np.loadtxt(rows,ndmin=2);d={k:data[:,i] for i,k in enumerate(names)}
    units=json.loads((OUT/'source-bindings.json').read_text())['constants']
    mass=d['star_mass']*units['G_SI']*units['Msun_kg']/1.3271244e20
    mass_ok=abs(mass/.197536385307-1)<=1e-8
    teff=10**d['log_Teff'];g=d['log_g']
    optical=mass_ok & (abs(teff-15800)<=300) & (abs(g-5.82)<=.15)
    score=((teff-15800)/100)**2+((g-5.82)/.05)**2
    eligible=np.flatnonzero(mass_ok);assert len(eligible)>0
    unconstrained_best=int(eligible[np.argmin(score[eligible])])
    candidates=np.flatnonzero(optical)
    best=int(candidates[np.argmin(score[candidates])]) if len(candidates) else unconstrained_best
    def row(i):
        return dict(model=int(d['model_number'][i]),star_age_yr=float(d['star_age'][i]),
          Teff_K=float(teff[i]),logg=float(g[i]),radius_source=float(10**d['log_R'][i]),
          mass_GM_sun=float(mass[i]),mass_relative_error=float(mass[i]/.197536385307-1),
          hydrogen_source_mass=float(d['total_mass_h1'][i]),period_days=float(d['period_days'][i]),
          diagnostic_distance_squared=float(score[i]),optical_candidate=bool(optical[i]))
    complete=(OUT/(name+'-execution.json')).exists()
    result=dict(classification='Counterexample candidate',run_completed=complete,
      last_model=int(d['model_number'][-1]),rows=len(data),matching_rows=int(optical.sum()),best=row(best),
      best_before_optical_cut=row(unconstrained_best),
      candidates=[row(int(i)) for i in np.flatnonzero(optical)],
      physical_GR_matching_completed=False,formation_history_and_mass_loss_rate_robustness=False)
    if complete:
        log=(folder/'execution.log').read_text()
        result['normal_termination']='termination code: max_model_number' in log or 'termination code: max_age' in log
        assert 'finished doing relax mass' in log
        index=np.loadtxt(folder/'LOGS1/profiles.index',skiprows=1,dtype=int)
        at=index[index[:,0]==result['best']['model']];assert len(at)==1
        path=folder/'LOGS1'/f'profile{at[0,2]}.data'
        dest=OUT/'mass-trial-best.data.gz'
        with dest.open('wb') as f:
            with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0) as g:
                with path.open('rb') as p:shutil.copyfileobj(p,g)
        result['best_raw_profile_sha256']=sha(path);result['best_profile']=str(path)
    save('mass-selection.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='candidates'},indent=2))


def retain():
    import numpy as np
    assert json.loads((OUT/'mass-selection.json').read_text())['run_completed']
    for rel in ['star/private/relax.f90','star/job/run_star_support.f90',
                'data/net_data/nets/basic.net','data/net_data/nets/add_hot_cno','data/net_data/nets/add_cno_extras']:
        p=OUT/'sources/mesa-r7624'/rel;p.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(CACHE/'mesa-r7624'/rel,p)
    def gzcopy(src,dest):
        dest.parent.mkdir(parents=True,exist_ok=True)
        with dest.open('wb') as f:
            with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0) as g:
                with src.open('rb') as p:shutil.copyfileobj(p,g)
    for name in ['restore19000','calcium19000','state19000','replay18000','mass18000']:
        record=json.loads((OUT/(name+'-inputs.json')).read_text());folder=CACHE/name
        dest=OUT/'runs'/name;dest.mkdir(parents=True,exist_ok=True)
        for rel in ['inlist1','cno_extras.net']:
            if rel not in record['sha256']:continue
            assert sha(folder/rel)==record['sha256'][rel]
            shutil.copy2(folder/rel,dest/rel)
        (dest/'restart-input.txt').write_text(record['photo']+'\n')
        assert sha(dest/'restart-input.txt')==record['sha256']['.restart']
        if (folder/'LOGS1/history.data').exists():gzcopy(folder/'LOGS1/history.data',dest/'history.data.gz')
        if (folder/'LOGS1/profiles.index').exists():shutil.copy2(folder/'LOGS1/profiles.index',dest/'profiles.index')
    gzcopy(CACHE/'state19000/LOGS1/profile1.data',OUT/'restored-model19000.data.gz')
    index=np.loadtxt(CACHE/'replay18000/LOGS1/profiles.index',skiprows=1,dtype=int)
    at=index[index[:,0]==19000];assert len(at)==1
    gzcopy(CACHE/'replay18000/LOGS1'/f'profile{at[0,2]}.data',OUT/'replayed-model19000.data.gz')
    src=PREV/'extracted/rotation_diffusion_3.4/LOGS1/history.data'
    expected=json.loads((ROOT/'outputs/thermal-wd18/profile-selection.json').read_text())['history_sha256']
    assert sha(src)==expected
    with (OUT/'published-history-window.data.gz').open('wb') as f:
        with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0) as g,src.open('rb') as original:
            for _ in range(6):g.write(next(original))
            for line in original:
                if line.strip() and 18000<=int(line.split(maxsplit=1)[0])<=19001:g.write(line)
    print('Retained source-derived reference window, run histories, actual controls and endpoint structures')


def maintain():
    trial=json.loads((OUT/'mass-selection.json').read_text());b=trial['best']
    assert trial['run_completed'] and trial['normal_termination']
    result=json.loads((OUT/'mass-response.json').read_text())
    mass_text=(f"사전 질량·광학 탐색 기준을 함께 만족한 행은 {trial['matching_rows']}개다. "
      f"기록한 최적 모델 {b['model']}은 Teff={b['Teff_K']:.6f} K, logg={b['logg']:.8f}, "
      f"GM/기준 GM_sun={b['mass_GM_sun']:.14f}이다.")
    additions={
      'model-definition':"분류: Proven. 공개 MESA 7624 사진에 필요한 ca40 핵종과 실행 중 유지되던 확산·대기 경계 설정을 명시적으로 복원했다. 원래 구조의 3,095개 셀·29개 구조 및 조성 열, 다음 한 단계, 18000→19000의 1,000단계 주요 이력과 최종 내부 구조가 일치했다. 원래 배포 상수도 직접 확인했다.\n\n분류: Counterexample candidate. 실제 외피 제거와 후속 열 진화로 Newtonian GM 목표에 맞춘 구조를 계산했다. "+mass_text+" 이는 GR 중력질량 및 형성 이력까지 검증한 matching이 아니다.",
      'observable-targets':"분류: Proven. 원래 저장되지 않았던 광학 후보 모델 18969의 실제 내부 구조를 재실행으로 확보했다. Teff=15786.146896 K, logg=5.82750698이다. 원래 단위의 GM 환산 질량은 기준보다 0.254509% 높다. 질량을 보정한 별도 진화의 탐색 결과는 다음과 같다. "+mass_text+" 관측 우도 적합이나 독립 반지름 측정으로 취급하지 않는다.",
      'adiabatic-limit':f"분류: Proven. 재현한 광학 모델 18969에서 정의한 고정 밀도 구각 scalar 모형은 eta=0.000165267742<1이고, 궤도 주파수의 상대 정적 차이 상계는 2.03086945e-10이다. 질량 보정 진화의 기록된 최적 구조에 동일한 조건부 식을 적용한 상계는 {result['relative_bound_upper']:.9e}이다. 두 결과 모두 실제 시간 의존 열 별의 전체 미분 오차 보장이 아니다.",
      'nonadiabatic-regime':"분류: Proven. 광학 후보 구조를 보간 없이 복원해도 영 배경·고정 밀도·평탄 시공간의 독립 scalar 모형은 궤도 주파수에서 정적 응답과 극히 가깝다. 실제 진화와 내부 구조의 수치 재현을 확인했지만 유체·metric 결합이나 중성자별 모드의 비단열 응답을 배제하지 않았다.",
      'failure-ledger-dynamic-chi':"분류: Proven. 초기 재시작 실패는 기본 cno_extras의 21핵종과 사진의 22핵종 불일치였다. ca40을 추가한 뒤에도 대기·확산 제어의 재시작 누락으로 다음 단계 반지름이 0.343% 달랐다. 이 설정을 원래 사용자 코드에 맞게 복원한 후 상태·한 단계·1,000단계 대조가 통과했다.\n\n분류: Counterexample candidate. "+mass_text+" 한 외피 제거율과 한 초기 이력만 검사했다. GR 질량 정의, 이력·제거율 및 격자 의존성, 실제 별의 전구간 오차 보장, 비영 배경의 전체 force/readout과 비선형 관측 추론은 별도 미완료다."}
    old=json.loads((ROOT/'outputs/thermal-wd18/manifest.json').read_text())['sha256']
    rels=['docs/'+x+'.md' for x in additions]+['paper/revision-manifest.json']
    dest=OUT/'request18-notes';dest.mkdir(exist_ok=False);bindings={}
    for rel in rels:
        assert sha(ROOT/rel)==old[rel]
        snap=dest/Path(rel).name;shutil.copy2(ROOT/rel,snap)
        bindings[rel]=dict(snapshot=snap.relative_to(ROOT).as_posix(),sha256=old[rel],historical_manifest=rel.startswith('paper/'))
    save('historical-note-bindings.json',bindings)
    revision=json.loads((ROOT/'paper/revision-manifest.json').read_text())
    for stem,body in additions.items():
        rel='docs/'+stem+'.md'
        with (ROOT/rel).open('ab') as f:f.write(('\n\n## Request 19 열 진화 재현과 질량 보정 대조\n\n'+body+'\n\n세부 근거: [한글 보고서](../notes/REQUEST19_THERMAL_RESTART_KO.md).\n').encode())
        revision['sha256'][rel]=sha(ROOT/rel)
    revision['request19_supporting_note_update']=dict(evidence_manifest='outputs/thermal-restart19/manifest.json',
      historical_notes='outputs/thermal-restart19/historical-note-bindings.json',
      status='MESA 진화 재현·광학 구조 복원·외피 질량 보정 대조 완료; GR 전체 matching·인증·추론 미완료',
      artifact_status='원고 PDF와 ZIP은 Request 12 동결본이다.')
    (ROOT/'paper/revision-manifest.json').write_text(json.dumps(revision,ensure_ascii=False,indent=2)+'\n')
    with (ROOT/'notes/REQUEST19_THERMAL_RESTART_KO.md').open('a') as f:
        f.write('\n## 외피 제거 후 진화 결과\n\n분류: Counterexample candidate. '+mass_text+
          f" 총 {trial['rows']}단계를 실행했다. 최대 제거율은 1e-6 원래 태양질량/년이다. "
          "MESA 내부 완화는 일시적으로 수치 제어를 바꾸고 내부 완화 뒤 나이·모델 번호를 복원하므로, 이 조작을 그대로 실제 쌍성의 형성 이력으로 해석하지 않는다.\n\n"+
          f"분류: Proven. 기록한 최적 구조의 고정 구각 scalar 상계는 {result['relative_bound_upper']:.9e}, "
          f"정적 응답은 {float(result['susceptibility_m']):.9f} m이다. 기호 전파식과 구간 연산을 재검증했다.\n\n"+
          "분류: Conjectural. 제거율·초기 외피 이력·격자 의존성, GR 질량 정의와 유체 결합, 실제 별의 전구간 미분 오차 보장, 전체 물리 timing과 비선형 관측 추론은 미완료다. 현재 묶음은 열 구조 수치 실험과 조건부 정리의 진전이며 제출 원고의 완성 판정은 아니다.\n")


def check():
    import numpy as np,mpmath as mp
    from thermal_wd import mesa,shell_response
    def hist(path):
        with gzip.open(path,'rt') as f:lines=f.readlines()
        data=np.loadtxt(lines[6:],ndmin=2)
        return {k:data[:,i] for i,k in enumerate(lines[5].split())}
    ref=hist(OUT/'published-history-window.data.gz');actual=hist(OUT/'runs/replay18000/history.data.gz')
    assert np.array_equal(actual['model_number'],np.arange(18001,19001))
    for k in ['star_mass','log_R','log_Teff','log_center_T','log_center_Rho','log_dt']:
        assert np.array_equal(ref[k][1:1001],actual[k]),k
    _,a=mesa(ROOT/'outputs/thermal-wd18/sources/profile388.data.gz')
    keys=['logT','logRho','logP','logR','mass','radius','omega']+json.loads((OUT/'network-diagnosis.json').read_text())['profile_isotopes']
    for filename in ['restored-model19000.data.gz','replayed-model19000.data.gz']:
        _,b=mesa(OUT/filename)
        assert all(np.array_equal(a[k],b[k]) for k in keys)
    trial=json.loads((OUT/'mass-selection.json').read_text());h,d=mesa(OUT/'mass-trial-best.data.gz')
    assert trial['run_completed'] and trial['normal_termination']
    assert int(h['model_number'])==trial['best']['model']
    assert abs(float(h['Teff'])/trial['best']['Teff_K']-1)<1e-14
    assert abs(float(h['photosphere_r'])/trial['best']['radius_source']-1)<1e-14
    mp.mp.dps=80
    for stem in ['target','mass']:
        data=np.load(OUT/(stem+'-shells.npz'));result=json.loads((OUT/(stem+'-response.json')).read_text())
        value,eta=shell_response(data['r_m'],data['m_geom_m'])
        lo,hi=map(mp.mpf,result['susceptibility_interval_m'].strip('[]').split(','))
        assert lo<=value<=hi and eta<1
        lo,hi=map(mp.mpf,result['relative_bound_interval'].strip('[]').split(','))
        assert mp.mpf(result['relative_bound_upper'])>=hi
        assert not result['full_stellar_or_gr_error_certificate']
    # Independent closed-form uniform sphere checks the shared propagator.
    v,_=shell_response(np.array([0.,1.]),np.array([0.,.001]))
    z=mp.sqrt(12*mp.mpf(.001));assert abs(v-(mp.tan(z)/z-1))<mp.mpf('1e-70')
    save('audit.json',dict(classification='Proven',frozen_history_rows_reproduced=1000,
      identical_endpoint_cells=3095,identical_structure_and_composition_columns=len(keys),
      selected_mass_profile_matches_history=True,independent_uniform_sphere=True,
      scalar_static_enclosures_recomputed_at_80_digits=True,outward_frequency_bounds_checked=True,
      full_stellar_error_certificate=False,complete_nonlinear_observational_inference=False))
    print('PASS: frozen 1000-row trajectory, both 3095-cell endpoints, mass-profile history, uniform sphere and scalar enclosures')


def seal():
    import numpy as np,mpmath,sympy,thermal_wd,stellar_matching,nonzero_drive
    trial=json.loads((OUT/'mass-selection.json').read_text())
    save('provenance.json',dict(before_task_checkpoint='15b77ea',before_mass_trial_checkpoint='4c630db',
      interpreter=sys.executable,versions={m.__name__:m.__version__ for m in [np,mpmath,sympy]},
      module_file_sha256={str(m.__file__):sha(m.__file__) for m in [np,mpmath,sympy,thermal_wd,stellar_matching,nonzero_drive]},
      input_sha256={'outputs/validated-variational/ivp.hex':sha(ROOT/'outputs/validated-variational/ivp.hex')},
      producer='verification/thermal_restart.py; see Korean report for ordered commands',
      archival_runtime_binding='environment.json contains verified full source/SDK archives and data/library SHA256'))
    save('gates.json',dict(classification='Proven',theorem_progress=True,
      original_thermal_evolution_reproduced=True,optical_epoch_structure_obtained=True,
      explicit_network_and_runtime_control_reconstruction_validated=True,
      envelope_mass_relaxation_and_followup_evolution_run=True,
      Newtonian_mass_and_optical_search_candidate_found=trial['matching_rows']>0,
      fixed_shell_scalar_response_certified=True,exact_J0337_GR_mass_and_structure_matched=False,
      original_hidden_network_file_acquired=False,full_stellar_interval_certificate=False,
      formation_history_and_mass_loss_rate_robustness=False,
      genuine_orbital_timescale_state_established=False,full_physical_dynamic_force_and_readout=False,
      full_span_variational_certificate=False,full_28_parameter_initialization_certificate=False,
      complete_nonlinear_observational_inference=False))
    paths=[p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    paths += [ROOT/'verification/thermal_restart.py',ROOT/'notes/REQUEST19_THERMAL_RESTART_KO.md']
    paths += [ROOT/x for x in json.loads((OUT/'historical-note-bindings.json').read_text())]
    save('manifest.json',dict(classification='Proven',scope='Request 19 열 진화 재현·외피 질량 보정·조건부 scalar 경계',
      sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(paths)}))


def verify():
    histories={
      'outputs/validated-variational/manifest.json':'outputs/remaining-levers15/historical-note-bindings.json',
      'outputs/remaining-levers15/manifest.json':'outputs/nbody-readout16/historical-note-bindings.json',
      'outputs/nbody-readout16/manifest.json':'outputs/nonzero-drive17/historical-note-bindings.json',
      'outputs/nonzero-drive17/manifest.json':'outputs/thermal-wd18/historical-note-bindings.json',
      'outputs/thermal-wd18/manifest.json':'outputs/thermal-restart19/historical-note-bindings.json'}
    count=0
    for label in ['outputs/research-remediation/manifest.json',*histories,'paper/revision-manifest.json','outputs/thermal-restart19/manifest.json']:
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
    print('PASS:',count,'현재·역사적 SHA, 원고 동결 및 전체 미완료 판정 보존')


if __name__=='__main__':globals()[sys.argv[1]]()
