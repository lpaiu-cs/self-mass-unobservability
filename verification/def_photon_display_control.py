"""Hold one calculation fixed; vary only its numerical-output energy window."""
from pathlib import Path
from http.cookiejar import CookieJar
import json
import signal
import time
import urllib.parse
import urllib.request
from bs4 import BeautifulSoup
import def_photon_edge_repair as previous

ex=previous.ex;h=previous.h;native=previous.native
OUT=previous.OUT.parent/'def-photon-display-control'


def main():
    assert not OUT.exists();OUT.mkdir();request=json.loads((previous.OUT/'requests.json').read_text())[0]
    fields=dict(request['fields'],mixname='phDisplay');ex.write(OUT/'input.json',fields)
    ex.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='0f6b09ed',
        claim='Distinguish generation loss from output-window filtering using one fresh physical query and two output requests in the same cookie session.',
        intervention='Keep the submit request fixed. Reuse its returned form unchanged once, then change only numerical-output egplow/egphigh outward by 1e-4 relative. Compare row count, overlapping opacity values and gray means; do not infer a physical recalculation from a display change.',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in [Path(__file__),OUT/'input.json']},budget=dict(submit=1,results=2,hard_seconds=20)))
    start=time.monotonic();signal.alarm(20);opener=urllib.request.build_opener(urllib.request.HTTPCookieProcessor(CookieJar()))
    def post(url,value,label):
        with opener.open(url,urllib.parse.urlencode(value).encode(),timeout=15) as response:raw=response.read()
        (OUT/(label+'.html')).write_bytes(raw);return BeautifulSoup(raw,'html.parser')
    try:
        soup=post('https://aphysics2.lanl.gov/submit',fields,'submit')
        form=next(f for f in soup.find_all('form') if f.find('input',attrs={'name':'output','value':'tabcol'}))
        second=native.retrieval.fields(form);results=[]
        for label,params in [('original',second),('outward',dict(second,egplow=f"{float(second['egplow'])*(1-1e-4):.17g}",egphigh=f"{float(second['egphigh'])*(1+1e-4):.17g}"))]:
            ex.write(OUT/(label+'-request.json'),params);page=post('https://aphysics2.lanl.gov/results',params,label)
            text=page.find('code').get_text('\n',strip=True).replace('\xa0',' ');(OUT/(label+'.txt')).write_text(text)
            lines=[s.strip() for s in text.splitlines() if s.strip()];n=int(next(s for s in lines if s.startswith('Photon grid')).split()[-2])
            j=next(i for i,s in enumerate(lines) if s.startswith('Energy') and 'density =' in s)
            results.append(dict(label=label,groups=n,rows=[s.split() for s in lines[j+1:j+1+n]]))
        ex.write(OUT/'result.json',dict(classification='Counterexample candidate',checks=results,seconds=time.monotonic()-start))
        print('DISPLAY CONTROL',results,flush=True)
    finally:signal.alarm(0)


if __name__=='__main__':main()
