"""Check the delivered submission and catch packaging the wrong checkout."""
import hashlib
import json
import re
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    manifest = json.loads((ROOT/'paper/submission-manifest.json').read_text())
    master = json.loads((ROOT/'paper/revision-manifest.json').read_text())
    for path, expected in manifest['sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == expected, path
        assert master['sha256'][path] == expected, path
    assert hashlib.sha256((ROOT/'paper/submission-manifest.json').read_bytes()).hexdigest() == master['sha256']['paper/submission-manifest.json']
    with zipfile.ZipFile(ROOT/'output/submission/free-fall-identifiability-source.zip') as archive:
        for member, source in [('main.tex', 'paper/main.tex'), ('references.bib', 'paper/references.bib'),
                               ('figures/conditional-intervals.pdf', 'paper/figures/conditional-intervals.pdf')]:
            assert archive.read(member) == (ROOT/source).read_bytes(), source
        bibkeys = set(re.findall(r'\\bibitem(?:\[[\s\S]*?\])?\{([^}]+)\}', archive.read('main.bbl').decode()))
        for source in ['manuscript.md', 'supplement.md']:
            text = (ROOT/'paper'/source).read_text(encoding='utf-8')
            cited = {key.strip() for group in re.findall(r'\\(?:cite|nocite)\{([^}]+)\}', text) for key in group.split(',')}
            assert cited <= bibkeys, sorted(cited-bibkeys)
        assert 'supplemental_material' in bibkeys
    for name, delivered in [('main', 'free-fall-identifiability'), ('supplement', 'free-fall-identifiability-supplement')]:
        built = ROOT/f'paper/build/{name}.pdf'
        if built.exists():
            assert built.read_bytes() == (ROOT/f'output/pdf/{delivered}.pdf').read_bytes(), name
    for path in ['paper/main.tex', 'paper/supplement.tex', 'output/submission/cover-letter-prd.tex']:
        text = (ROOT/path).read_text(encoding='utf-8')
        assert 'lpaiu.cs@gmail.com' in text and 'Independent researcher' in text, path
        assert '\\fillin' not in text, path
    print(f'PASS: {len(manifest["sha256"])} submission bindings; ZIP/source identity; complete main/SM references; author metadata')


if __name__ == '__main__':
    main()
