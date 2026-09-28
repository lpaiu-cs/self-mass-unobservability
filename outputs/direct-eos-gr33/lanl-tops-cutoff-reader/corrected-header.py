def header(path,request):
    lines=[x.strip() for x in path.read_text().splitlines() if x.strip()]
    fields=request['fields'];tokens=fields['mixture'].split()
    amounts={int(tokens[i+1]):(F(Decimal(tokens[i])),F(Decimal(tokens[i+2]))) for i in range(0,len(tokens),3)}
    nt=sum(x[0] for x in amounts.values());mt=sum(x[0]*x[1] for x in amounts.values())
    start=next(i for i,x in enumerate(lines) if x.startswith('No. Fraction'))+1;echo=[];ids={}
    while not lines[start].startswith('Temperature grid'):
        number,mass,Z,symbol,mat=lines[start].split();Z=int(Z);n,w=amounts[Z]
        echo.extend([audit.score(n/nt,number),audit.score(n*w/mt,mass)]);ids[Z]=int(mat);start+=1
    assert set(ids)==set(amounts)
    expected=next(x for x in json.loads((audit.OUT/'result.json').read_text())['rows'] if x.get('accepted'))['material_ids']
    assert ids=={x['Z']:x['returned'] for x in expected}
    T=lines[start+1].split();rho=lines[start+3].split();assert len(rho)==1
    Ts=sorted(F(Decimal(x)) for x in fields['temps'].split());assert len(Ts)==len(T)
    echo.extend(audit.score(a,b) for a,b in zip(Ts,sorted(T,key=Decimal)))
    density=F(Decimal(fields['dens']));echo.append(audit.score(density,rho[0]));means={}
    for i,line in enumerate(lines):
        if line.startswith('Density') and 'T=' in line:
            t=line.split('T=')[1].strip();v=lines[i+1].split();assert len(v)==5
            echo.append(audit.score(density,v[0]));means[t]=list(map(float,v))
    assert set(means)==set(T)
    groups=[];energies=[]
    if fields['datype']=='groups':
        i=next(i for i,x in enumerate(lines) if x.startswith('Photon grid'))
        count=int(re.search(r'(\d+) points',lines[i]).group(1));i+=1
        while not lines[i].startswith('Rosseland'):energies.extend(lines[i].split());i+=1
        assert len(energies)==count
        submitted=fields['energies'].split();assert len(submitted)==count+1
        second=json.loads((retrieval.OUT/(request['name']+'-results-request.json')).read_text())
        assert F(Decimal(second['egplow']))==F(Decimal(submitted[0])) and F(Decimal(second['egphigh']))==F(Decimal(submitted[-1]))
        echo.extend(audit.score(F(Decimal(a)),b) for a,b in zip(submitted[:-1],energies))
        i=next(i for i,x in enumerate(lines) if x.startswith('Energy') and 'density =' in x)
        tg,rg=lines[i].split('=')[1].split();echo.append(audit.score(density,rg))
        echo.append(audit.score(Ts[0],tg));groups=np.array([x.split() for x in lines[i+1:]],float)
        assert groups.shape==(count,3) and np.all(groups>0) and np.all(np.isfinite(groups))
        assert np.array_equal(groups[:,0],np.array(energies,float))
    passed=max(echo)<=1 and not any('warning' in x.lower() for x in lines)
    return dict(input_passed=passed,input_max_score=max(echo),means=means,groups=np.asarray(groups),
                mass_sum=float(mt),number_sum=float(nt))
