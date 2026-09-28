"""Package the focused paper and Supplemental Material from this checkout.

Copies the built PDFs, writes the plain-text abstract from paper/manuscript.md, and writes the journal source archive (main.tex,
main.bbl, references.bib, the one figure the focused paper uses, README.txt).
Refreshes submission hashes without dropping any existing research-input bindings.
"""
import hashlib, json, re, shutil, zipfile
from pathlib import Path
root = Path(__file__).resolve().parents[4]
documents = (('main', 'manuscript'), ('supplement', 'supplement'))
for name, source in documents:
    src = root/f'paper/build/{name}.pdf'
    inputs = [root/f'paper/{name}.tex', root/f'paper/{source}.md', root/'paper/references.bib']
    inputs += list((root/'paper/figures').glob('*.pdf'))
    assert src.stat().st_mtime > max(p.stat().st_mtime for p in inputs), f'rebuild {name}.pdf after its inputs'
cover = root/'output/submission/cover-letter-prd.pdf'
assert cover.stat().st_mtime > cover.with_suffix('.tex').stat().st_mtime, 'rebuild cover-letter-prd.pdf'
for name, dst in (('main', 'free-fall-identifiability.pdf'), ('supplement', 'free-fall-identifiability-supplement.pdf')):
    src = root/f'paper/build/{name}.pdf'
    shutil.copyfile(src, root/'output/pdf'/dst)

# Plain-text abstract for the submission form.
ms = (root/'paper/manuscript.md').read_text(encoding='utf-8')
abstract = ms.split('## Abstract', 1)[1].split('\n---', 1)[0].strip()
for a, b in ((r'\(K\)', 'K'), (r'\(N\)', 'N'), (r'\(N\ge2K-1\)', 'N ≥ 2K − 1'), ('--', '–'), ('**', '')):
    abstract = abstract.replace(a, b)
assert not re.search(r'[\\$*]', abstract), abstract
(root/'output/submission/abstract-plain.txt').write_text(abstract + '\n', encoding='utf-8', newline='\n')

title = ms.split('\n', 1)[0][2:].strip()
readme = (f'{title}\n'
          'Juneyoung Kim, manuscript of 28 September 2026.\n\n'
          'Compile main.tex with latexmk -pdf main.tex (pdfLaTeX and BibTeX) or tectonic main.tex.\n'
          'The bibliography (references.bib), the generated main.bbl and the figure are included.\n'
          'The Supplemental Material is submitted separately as a PDF (free-fall-identifiability-supplement.pdf).\n'
          'main.tex is generated from manuscript.md in https://github.com/lpaiu-cs/self-mass-unobservability (paper/build_manuscript.py);\n'
          'data, code and SHA-256 manifests are in that repository, as stated in the Data and code availability section.\n')
members = {'paper/main.tex': 'main.tex', 'paper/build/main.bbl': 'main.bbl', 'paper/references.bib': 'references.bib',
           'paper/figures/conditional-intervals.pdf': 'figures/conditional-intervals.pdf'}
tex = (root/'paper/main.tex').read_text(encoding='utf-8')
assert re.findall(r'\\includegraphics\[[^]]*\]\{([^}]*)\}', tex) == ['figures/conditional-intervals.pdf']
archive = root/'output/submission/free-fall-identifiability-source.zip'
with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED) as z:
    for src, dst in members.items(): z.write(root/src, dst)
    z.writestr('README.txt', readme)
print('abstract words', len(abstract.split()), 'zip', archive.stat().st_size, sorted(zipfile.ZipFile(archive).namelist()))

# Bind the delivered files while preserving all earlier research inputs and receipts.
paths = list(members) + [
    'paper/manuscript.md', 'paper/supplement.md', 'paper/supplement.tex', 'paper/build_manuscript.py',
    'paper/README.md', 'README.md', 'CITATION.cff', 'PUBLIC_SNAPSHOT.md',
    'output/pdf/free-fall-identifiability.pdf', 'output/pdf/free-fall-identifiability-supplement.pdf',
    'output/submission/cover-letter-prd.tex', 'output/submission/cover-letter-prd.pdf',
    'output/submission/abstract-plain.txt', 'output/submission/free-fall-identifiability-source.zip',
    'output/submission/submission-checklist-prd.md', 'verification/replay_public_inference.py',
    'verification/check_submission_package.py', 'outputs/research-completion/public-inference.json',
    'notes/REQUEST294_PRE_SUBMISSION_AUDIT_KO.md', 'notes/REQUEST295_SUBMISSION_CORRECTIONS_KO.md',
    'docs/model-definition.md', 'docs/observable-targets.md', 'docs/adiabatic-limit.md',
    'docs/nonadiabatic-regime.md', 'docs/failure-ledger-dynamic-chi.md',
    Path(__file__).resolve().relative_to(root).as_posix()]
# main.bbl is delivered inside the ZIP; paper/build is an ignored local build directory.
paths.remove('paper/build/main.bbl')
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
submission = dict(revision='prd-submission-2026-09-28', reviewed_commit='6a0279b1276d60d9440f12e211a1f7e38daf7d05',
                  classification='Imported from prior work',
                  scope='Submission corrections and public replay; no new empirical SEP bound or full-goal completion.',
                  sha256={p: sha(root/p) for p in paths})
manifest = root/'paper/submission-manifest.json'
manifest.write_text(json.dumps(submission, indent=2) + '\n', encoding='utf-8', newline='\n')
master = root/'paper/revision-manifest.json'
text = master.read_text(encoding='utf-8'); data = json.loads(text)
data['sha256'].update(submission['sha256'])
data['sha256']['paper/submission-manifest.json'] = sha(manifest)
data['submission_revision'] = {k: v for k, v in submission.items() if k != 'sha256'}
master.write_text(json.dumps(data, indent=2 if text.startswith('{\n  "') else 1, ensure_ascii=False) + '\n',
                  encoding='utf-8', newline='\n')
print('submission manifest', len(submission['sha256']), 'bindings; preserved master', len(data['sha256']), 'bindings')
