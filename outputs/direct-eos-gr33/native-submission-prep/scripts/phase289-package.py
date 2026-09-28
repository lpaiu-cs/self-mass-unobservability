"""Phase 289: package the submission outputs without touching paper/revision-manifest.json.

paper/package_revision.py rewrites revision-manifest.json from its fixed 9 September list and would drop thousands of later
bindings, so it is not used. This script copies the built PDF and writes a journal source archive (main.tex, main.bbl,
references.bib, figures, README.txt); the phase publisher binds the results.
"""
import shutil, zipfile
from pathlib import Path
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
pdf_src, pdf_dst = root/'paper/build/main.pdf', root/'output/pdf/free-fall-identifiability.pdf'
assert pdf_src.stat().st_mtime > (root/'paper/main.tex').stat().st_mtime, 'rebuild the PDF after main.tex'
shutil.copyfile(pdf_src, pdf_dst)
members = {'paper/main.tex': 'main.tex', 'paper/build/main.bbl': 'main.bbl', 'paper/references.bib': 'references.bib',
           **{f'paper/figures/{n}': f'figures/{n}' for n in ('conditional-intervals.pdf', 'coverage-validation.pdf',
                                                             'comparator-phase-validation.pdf', 'remaining-levers-validation.pdf')}}
readme = ('Static response and dynamical identifiability in free-fall tests: finite-order boundaries and a pulsar-triple application\n'
          'Juneyoung Kim, manuscript of 28 September 2026.\n\n'
          'Compile main.tex with latexmk -pdf main.tex (pdfLaTeX and BibTeX) or tectonic main.tex.\n'
          'The bibliography (references.bib), the generated main.bbl and the four figures are included.\n'
          'main.tex is generated from manuscript.md in https://github.com/lpaiu-cs/self-mass-unobservability (paper/build_manuscript.py);\n'
          'data, code and SHA-256 manifests are in that repository, as stated in the Data and code availability section.\n')
archive = root/'output/submission/free-fall-identifiability-source.zip'
with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED) as z:
    for src, dst in members.items(): z.write(root/src, dst)
    z.writestr('README.txt', readme)
print('pdf', pdf_dst.stat().st_size, 'zip', archive.stat().st_size, sorted(zipfile.ZipFile(archive).namelist()))
