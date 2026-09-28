eos=EOS();checks=[]
for xx,tt in [(-.01,13000),(-.11,14000),(-.37,9000),(-1.1,6900),(-2.3,3500),(-3.1,1900),(-4.2,850),(-5.3,390),(-5.9,260),(-3.7,12000),(-5.5,7000),(-.8,18000)]:
    t=np.log(tt);a=native(xx,t);p,u,gamma,T,kap,cv,s=eos.evaluate(np.array([np.exp(xx)]),np.array([t]))
    h=1e-4;rp=native(xx+h,t);rm=native(xx-h,t);tp=native(xx,t+h);tm=native(xx,t-h)
    ar=(rp-rm)/(2*h);at=(tp-tm)/(2*h);g=ar[1]/a[1]+at[1]/a[1]*(a[1]/a[0]-ar[2])/at[2]
    errors=[float(abs(p[0]*eos.rho0*C*C/a[1]-1)),float(abs(u[0]*C*C/a[2]-1)),float(abs(gamma[0]/g-1)),float(abs(cv[0]*C*C/at[2]-1))]
    laws=[float((tt*at[3]-at[2])/at[2]),float((tt*ar[3]-ar[2]+a[1]/a[0])/(a[1]/a[0]))]
    checks.append(dict(x=xx,T=tt,relative=errors,first_law=laws))
result=dict(classification='Counterexample candidate',passed=bool(max(max(c['relative']) for c in checks)<.002 and max(max(abs(np.array(c['first_law']))) for c in checks)<1e-4),
    checks=checks,EOS_calls=native.ion.calls,seconds=time.monotonic()-start,source_sha256=cold.sha(__file__))
write(OUT/'eos.json',result);print(json.dumps(result),flush=True);assert result['passed']
