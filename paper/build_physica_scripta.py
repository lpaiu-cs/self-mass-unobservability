"""Build the bounded Physica Scripta revision without modifying frozen editions.

Uses the existing Markdown builder, installed Pandoc and installed Tectonic.
No timing, stellar or Monte Carlo calculations are launched.
"""
import hashlib
import json
import os
import re
import shutil
import subprocess
import zipfile
from pathlib import Path

from build_manuscript import build_tex, parse_document

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/submission/physica-scripta'
OLD = ROOT / 'output/submission/grg'
BUILD = ROOT / 'tmp/pdfs/physica-scripta-build'
TECTONIC = os.environ.get('TECTONIC', 'tectonic')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compile_tex(path):
    subprocess.run([TECTONIC, '--keep-logs', '--keep-intermediates',
                    '--outdir', str(BUILD), str(path)], cwd=OUT, check=True)
    log = (BUILD / f'{path.stem}.log').read_text(errors='replace')
    assert not re.search(r'Undefined control sequence|undefined references|undefined citations|Overfull \\hbox', log), path


def main():
    frozen = json.loads((OLD / 'manifest.json').read_text())
    for name, expected in frozen['sha256'].items():
        assert digest(ROOT / name) == expected, f'frozen GRG file changed: {name}'
    BUILD.mkdir(parents=True, exist_ok=True)
    manuscript = (OUT / 'manuscript.md').read_text(encoding='utf-8')
    title, metadata, abstract, body = parse_document(manuscript)
    # Keep this journal's layout local; the older default builder stays unchanged.
    for source, stem in [('manuscript.md', 'main'), ('supplement.md', 'supplement')]:
        text = (OUT / source).read_text(encoding='utf-8')
        parsed = parse_document(text)
        tex = build_tex(*parsed, source_name=source).replace('[11pt]{article}', '[12pt]{article}')
        tex = tex.replace('% Edit the markdown source and rerun the build script.',
                          '% Source: output/submission/physica-scripta/' + source)
        tex = tex.replace('from paper/' + source, 'from output/submission/physica-scripta/' + source)
        # Supplemental historical tables were designed at 11pt; keep them legible
        # within 12pt text without shrinking the scientific prose or equations.
        if stem == 'supplement':
            tex = tex.replace('\\usepackage{graphicx,booktabs}',
                              '\\usepackage{graphicx,booktabs}\n\\usepackage{etoolbox}\n'
                              '\\AtBeginEnvironment{tabular}{\\small}')
        path = OUT / f'{stem}.tex'
        path.write_text(tex, encoding='utf-8')
        compile_tex(path)
        bibliography = BUILD / f'{stem}.bbl'
        shutil.copy2(bibliography, OUT / bibliography.name)
        if stem == 'main':
            # No external figures in the main article: inline the generated
            # bibliography so the native standalone editor can also compile it.
            tex = tex.replace('\\bibliographystyle{unsrtnat}\n\\bibliography{references}',
                              bibliography.read_text(encoding='utf-8'))
            path.write_text(tex, encoding='utf-8')
            compile_tex(path)
        shutil.copy2(BUILD / f'{stem}.pdf', OUT / ('manuscript.pdf' if stem == 'main' else 'supplement.pdf'))
    compile_tex(OUT / 'cover-letter.tex')
    shutil.copy2(BUILD / 'cover-letter.pdf', OUT / 'cover-letter.pdf')

    plain = re.sub(r'\*\*(?:Proven|Imported from prior work|Counterexample candidate|Conjectural)\.\*\*\s*', '', abstract)
    plain = plain.replace(r'\(', '').replace(r'\)', '').replace(r'\times', ' x ')
    plain = plain.replace(r'\times', ' x ').replace(r'\tau_\chi', 'tau_chi')
    plain = plain.replace(r'\ge', '>=').replace(r'\beta', 'beta')
    plain = plain.replace('\\', '').replace('{', '').replace('}', '')
    assert len(plain.split()) <= 300
    (OUT / 'abstract-plain.txt').write_text(plain + '\n', encoding='utf-8')
    fields = (f'Target: Physica Scripta\nArticle type: Research Paper\nTitle: {title}\n'
              'Author: Juneyoung Kim\nGiven name: Juneyoung\nFamily name: Kim\n'
              'Affiliation: Independent researcher\nCity: Seoul\nCountry: Republic of Korea\n'
              'Corresponding email: lpaiu.cs@gmail.com\n'
              'Keywords: dynamical identifiability; relaxation response; nuisance projection; pulsar timing; tests of gravity\n'
              'External funding: None\nCompeting interests: None\n'
              'Preferred route: Subscription publication, no APC\nPeer review: Single-anonymous\n'
              'Previous submission: GRG editorial rejection, 29 September 2026; no referee reports\n'
              'Supplement title: Supplementary Information\n'
              'Supplement description: Detailed derivations, archived timing audits and conditional stellar calculation; '
              'historical numbering retained, with unresolved derivative, noise and physical-matching limitations.\n'
              'Code title: Supplementary Code\n'
              'Code description: Exact projected-separation controls, frozen comparator records and compact inference replay; '
              'no timing engine, large runtime arrays or new astrophysical validation.\n'
              'Author policy declarations and final submission: to be completed personally by the submitting author.\n')
    (OUT / 'submission-fields.txt').write_text(fields, encoding='utf-8')

    with zipfile.ZipFile(OUT / 'manuscript-source.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
        for name in ['main.tex', 'main.bbl', 'references.bib']:
            archive.write(OUT / name, name)
    record = json.loads((ROOT / 'outputs/research-completion/public-inference.json').read_text())
    code_files = sorted(set(record['source_sha256']) | {
        'verification/projected_identifiability.py', 'verification/replay_public_inference.py',
        'symbolic/state_identifiability.py', 'outputs/research-completion/public-inference.json',
        'outputs/research-completion/comparator-audit.json',
    })
    for name, expected in record['source_sha256'].items():
        assert digest(ROOT / name) == expected, name
    code_readme = '''# Supplementary Code

Description: Exact projected-separation controls, frozen comparator records and compact inference replay; no timing engine, large runtime arrays or new astrophysical validation.

Companion: Finite-frequency identifiability of relaxation responses in free-fall timing (30 September 2026).

Requirements: Python 3.11+, NumPy and SymPy. Run from the extracted archive root:

    python verification/projected_identifiability.py
    python verification/replay_public_inference.py
    python symbolic/state_identifiability.py

The first command reports 15 exact controls, including coefficient distance 117/1000, erased information 0 and attenuated information 117*epsilon**2/1000. These are constructed dimensionless examples (Proven).

The second command reports the frozen threshold 12.8241766 (to the displayed precision), omnibus statistic about 16.3525 and six empty physical sections (Imported from prior work). It checks source hashes and rebuilds the physical columns. It is conditional on the exported coefficients, precision and threshold. It does not certify derivatives, repeat the noise calibration, fit a nonlinear timing model or produce an empirical SEP bound.

The third command reports 14 exact identities for positive moments, close poles, invisible filters and equilibrium/dynamics ambiguity. It rewrites only its derived analytic JSON in outputs/research-completion/.

outputs/research-completion/comparator-audit.json supplies every retained-information case behind main Table 2. For a comparator, take the minimum relative_information across its stored rows. Unit carrier amplitudes, three phase sets and three fixed covariance amplitudes are part of this archived calculation. Degree-five floating-point remnants are not positive information; the exact interpolation identity fixes zero.

Additional source files are included to satisfy the replay's provenance hash checks. They are not instructions to launch the omitted runtime or rerun the original calibration. Large binary arrays are not included. See the paper and Supplementary Information for remaining numerical and physical limitations.
'''
    with zipfile.ZipFile(OUT / 'supplementary-code.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
        archive.writestr('README.md', code_readme)
        for name in code_files:
            archive.write(ROOT / name, name)
        archive.writestr('manifest.json', json.dumps({name: digest(ROOT / name) for name in code_files}, indent=2) + '\n')
    # Reopen the delivered archives and check their actual bytes.
    for zip_name, paths in [('manuscript-source.zip', ['main.tex', 'main.bbl', 'references.bib']),
                            ('supplementary-code.zip', code_files)]:
        with zipfile.ZipFile(OUT / zip_name) as archive:
            for name in paths:
                source = OUT / name if zip_name == 'manuscript-source.zip' else ROOT / name
                assert archive.read(name) == source.read_bytes(), name
    files = [p for p in OUT.rglob('*') if p.is_file() and p.name != 'manifest.json']
    result = dict(
        revision='physica-scripta-2026-09-30', status='Prepared; not submitted',
        scope='Methods-focused revision; exact analytic controls and reused conditional numerical results; no empirical SEP bound.',
        abstract_word_count=len(plain.split()),
        source_sha256={str(p.relative_to(ROOT)).replace('\\', '/'): digest(p)
                       for p in [ROOT / 'paper/build_physica_scripta.py', ROOT / 'verification/projected_identifiability.py']},
        frozen_grg_manifest_sha256=digest(OLD / 'manifest.json'),
        sha256={str(p.relative_to(ROOT)).replace('\\', '/'): digest(p) for p in sorted(files)},
    )
    (OUT / 'manifest.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    for name, expected in result['sha256'].items():
        assert digest(ROOT / name) == expected, name
    print(f'PASS: {len(files)} package files, portable ZIP bytes, frozen GRG bindings; abstract {len(plain.split())} words')


if __name__ == '__main__':
    main()
