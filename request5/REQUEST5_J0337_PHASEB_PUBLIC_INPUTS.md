# Request 5 Phase B: Public Inputs Scout

`Request 5 Phase A` already translated published `J0337` SEP bounds into the
`(\sigma_1, \sigma_2)` EFT plane.

The next bounded question for **Phase B** is simpler:

> are the actual `TOA + code` inputs for a `J0337` three-body timing fit public
> and reproducibly acquirable?

This note closes that question.

## Artifacts

- scout script:
  [request5_j0337_phaseB_public_inputs.py](request5_j0337_phaseB_public_inputs.py)
- machine summary:
  [request5_j0337_phaseB_public_inputs_summary.json](request5_j0337_phaseB_public_inputs_summary.json)
- visual summary:
  [request5_j0337_phaseB_public_inputs_summary.svg](request5_j0337_phaseB_public_inputs_summary.svg)
- asset manifest:
  [request5_j0337_phaseB_public_inputs_manifest.tsv](request5_j0337_phaseB_public_inputs_manifest.tsv)

## Primary Sources

The public-input path is now anchored by primary sources, not just secondary
mentions.

### 2020 release

- paper / open-access description:
  [A&A 638, A24 (2020)](https://www.aanda.org/articles/aa/full_html/2020/06/aa38104-20/aa38104-20.html)
- Zenodo release:
  [10.5281/zenodo.3778978](https://zenodo.org/records/3778978)

Zenodo explicitly states that the release contains:

- source code of `Nutimo`,
- dataset,
- and results

for the 2020 `J0337` SEP paper.

### 2025 release

- paper / footnote pointing to data and code:
  [A&A 693, A143 (2025)](https://www.aanda.org/component/article?access=doi&doi=10.1051%2F0004-6361%2F202452100)
- Zenodo release:
  [10.5281/zenodo.13899771](https://zenodo.org/records/13899771)

This release explicitly says it contains:

- `Nutimo` source code,
- the dataset,
- analysis products,
- and results.

Guillaume Voisin's public page independently points to the same `Nutimo`
release family as the code path developed for `PSR J0337+1715`:
[My Work – Guillaume Voisin, Astrophysics](https://people-lux.obspm.fr/gvoisin/mywork/).

## Local Acquisition Result

The bounded scout mirrored the following locally:

- 2020:
  `README`, `nutimo.tar.bz2`, `Data_and_Results.tar.bz2`
- 2025:
  `Readme.md`, `nutimo.tar.bz2`, `Data.tar.bz2`

The very large 2025 `Analysis.tar.bz2` remains remote-only in this bounded
scout because it is not required to answer the public-input question.

## What Was Found

### Public TOA files are present

Inside the mirrored bundles:

- 2020 data bundle contains timing-data products such as
  `f4t1200s-FDM-20190514-gt600s-clean.tim-sorted-supershort`,
  `4ft1200s-20191018-FDM-gt600s-superclean-no_sun_5deg.tim`,
  and `0337+1715.par-nofit-DE430`
- 2025 data bundle contains
  `0337_20211005-sorted-sun5deg-res25microsec.tim`
  and
  `0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim`

The 2025 `Data/Readme.md` says explicitly that the clipped file is the one used
in the paper.

So `Phase B` is no longer blocked on "no public TOA input".

### Public three-body timing code is present

Both releases include `nutimo.tar.bz2`.

The 2025 code release contains:

- `nutimo/src/AllTheories3Bodies.cpp`
- `nutimo/src/Fittriple-init.cpp`
- `nutimo/src/Fittriple-compute.cpp`
- `nutimo/src/MCMC_parallelaffineinvariant.cpp`
- `nutimo/src/makefile-original`
- `nutimo/install_script.sh`
- Cython wrappers and Python analysis helpers

So `Phase B` is also no longer blocked on "no public timing model code".

## Current Boundary

The new bottleneck is not data access.
It is local build/runtime closure.

The 2025 `Nutimo` README and install script make the dependency stack explicit:

- `Tempo2`
- `Minuit`
- `Boost`
- `Acor`
- `python`
- and Cython extensions

The current host probe shows:

- `g++`, `gcc`, `python3` exist
- but `tempo2`, `root-config`, and `cython` do not

So the most accurate current statement is:

> `Request 5 Phase B` now has a real public-input path, but not yet a closed
> local runtime path.

That runtime question has now also been checked in a bounded way:
`REQUEST5_J0337_PHASEB_BUILD_FEASIBILITY.md` shows that the current host fails
first on Apple clang `-fopenmp`, and then on missing `boost/numeric/odeint.hpp`
once `g++-14` is forced. So the local blocker is explicit dependency/runtime
closure, not uncertainty about public inputs.

## Project Consequence

This changes the status of `Request 5` materially:

- `Phase B` is no longer hypothetical in terms of data/code availability,
- the next gate is a bounded `Nutimo` build/runtime feasibility test,
- and if that gate fails, the project can still say exactly where the public
  path stops.

That is a much stronger position than the earlier "TOA and code may or may not
be public" state.

## Status Line

The most accurate one-line status is now:

> `J0337 Phase B` has a real public-input branch: public `TOA`, `par`, and
> `Nutimo` code are locally mirrored, so the next bounded question is only
> whether the current workspace can close the required `Tempo2/Minuit/Boost`
> runtime stack without turning into an infrastructure project.
