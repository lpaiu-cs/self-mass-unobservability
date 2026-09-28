"""Final pre-submission checks: bibliography (DOI/arXiv resolution, cited keys, printed notes), spelling variants, plain 'pi'/'c_chi',
abstract length, unresolved references in the log."""
import re, urllib.request, json
from pathlib import Path
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
bib = (root/'paper/references.bib').read_text(encoding='utf-8')
md = (root/'paper/manuscript.md').read_text(encoding='utf-8')
bbl = (root/'paper/build/main.bbl').read_text(encoding='utf-8')
entries = re.findall(r'@(\w+)\{([^,]+),(.*?)\n\}', bib, re.S)
cited = sorted(set(k.strip() for c in re.findall(r'\\cite\{([^}]+)\}', md) for k in c.split(',')))
print('bib entries', len(entries), 'cited keys', len(cited), 'uncited', sorted(set(k for _, k, _ in entries) - set(cited)), 'missing', sorted(set(cited) - set(k for _, k, _ in entries)))
def head(url):
    req = urllib.request.Request(url, method='HEAD', headers={'User-Agent': 'Mozilla/5.0 (submission check)'})
    try:
        with urllib.request.urlopen(req, timeout=20) as r: return r.status
    except Exception as e: return getattr(e, 'code', str(e)[:40])
for typ, key, body in entries:
    doi = re.search(r'doi\s*=\s*\{([^}]+)\}', body); ep = re.search(r'eprint\s*=\s*\{([^}]+)\}', body); url = re.search(r'url\s*=\s*\{([^}]+)\}', body)
    note = re.search(r'note\s*=\s*\{(.*?)\}\s*,?\s*$', body, re.S)
    res = []
    if doi:
        try:
            with urllib.request.urlopen(f'https://api.crossref.org/works/{doi.group(1)}', timeout=20) as r:
                t = json.load(r)['message']['title'][0]; res.append('crossref:' + t[:60])
        except Exception as e: res.append('crossref:FAIL ' + str(getattr(e, 'code', e))[:30])
    if ep: res.append('arXiv:' + str(head('https://arxiv.org/abs/' + ep.group(1))))
    if url and not doi: res.append('url:' + str(head(url.group(1))))
    if not (doi or ep or url): res.append('NO IDENTIFIER')
    print(f'{key:30s}', ' | '.join(res), '| note' if note else '')
for word in ['Source of the', 'draft', 'TODO', 'REQUEST', 'repository ledger']:
    if word in bbl: print('bbl contains:', word)
for pat in [r'\bcentre\b', r'\bmetres?\b', r'behaviour', r'modell', r'analys(e|ed|ing)\b', r'normalis', r'colour', r'\bc_chi\b', r'\bpi\b', r'Nancay', r'\d-\d{4}\b', r'\d\.\d+-\d']:
    hits = [(i + 1, l[max(0, m.start() - 30):m.end() + 20]) for i, l in enumerate(md.split('\n')) for m in re.finditer(pat, l)]
    if hits: print(pat, len(hits), hits[:6])
abstract = re.search(r'## Abstract\s*\n\s*\n(.*?)\n\s*\n', md, re.S).group(1)
plain = re.sub(r'\\\(|\\\)|\*\*', '', abstract); print('abstract words', len(plain.split()), 'chars', len(plain))
log = (root/'paper/build/main.log').read_text(encoding='utf-8', errors='replace')
print('log undefined:', log.count('undefined'), 'overfull:', log.count('Overfull'), 'underfull:', log.count('Underfull'))
