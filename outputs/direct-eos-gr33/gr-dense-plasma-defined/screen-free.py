def free_mp(rs,ge,z):
    # A decimal real literal without d0 is rounded as Fortran default real.
    q=lambda s:mp.mpf(float(np.float32(s)))
    X=q('.0140047')/rs;logz=mp.log(z);z13=mp.exp(logz/3)
    cdh=z/q('1.73205')*(mp.sqrt(z+1)**3-mp.sqrt(z)**3-1)
    ctf=z*z*q('.2513')*(z13-1+q('.2')/mp.sqrt(z13))
    p01=q('1.11')*mp.exp(q('.475')*logz);p03=q('.2')+q('.078')*logz**2;power=q('1.16')+q('.08')*logz
    tx=ge**power
    cor1=1+(z-1)/9*rs**3/(1+6*rs**2)*(1+1/(mp.mpf(float('.001'))*z*z+2*ge))
    cor0=1+q('.78')*mp.sqrt(ge/z)*rs**3/(ge*z**3+21*rs**3)
    h1=(1+X*X/5)/(1+q('.18')/mp.sqrt(mp.sqrt(z))*X+(q('.2')+q('.37')/mp.sqrt(z))*X*X)
    up=cdh*mp.sqrt(ge)+p01*ctf*tx*cor0*h1
    den=1+(p03*mp.sqrt(ge)+p01/rs*tx*cor1)/mp.sqrt(1+X*X)
    return -up/den*ge
