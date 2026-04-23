# Request 4: External LLR Estimator Scout

`Request 4` has now crossed its `MLRS` stop rule:

- bundled two-case replay works,
- bounded weak-field observation-operator seams work,
- public monthly `CRD` promotion reaches parser/recalc,
- but the public `raw frd -> usable frcal / 93 stream` layer was not found.

So the final weak-field posterior must now be handed off to a more mature
external estimator path.

This note ranks the best candidates found in the bounded scout.

## Artifacts

- scout script:
  [request4_llr_external_estimator_scout.py](request4_llr_external_estimator_scout.py)
- machine summary:
  [request4_llr_external_estimator_summary.json](request4_llr_external_estimator_summary.json)
- visual summary:
  [request4_llr_external_estimator_summary.svg](request4_llr_external_estimator_summary.svg)
- candidate table:
  [request4_llr_external_estimator_candidates.tsv](request4_llr_external_estimator_candidates.tsv)

## Ranked Outcome

### Tier A: `PEP`

`PEP` is the strongest external hand-off candidate found.

Why:

- the `PEP` literature describes it as an open-source general-purpose
  astrometric data-analysis program,
- the same literature explicitly states that source code, utilities,
  documentation, and unit tests are publicly available via `GitLab`,
- `ASCL` lists `PEP` as a code used for solar-system astrometric data analysis,
- and `ILRS` literature shows it being used on real `LLR` analysis problems.

This makes `PEP` the best match to the actual gap left by the frozen `MLRS`
branch:

- not just an ephemeris reader,
- not just a data archive,
- but a real fit engine with `LLR` pedigree.

Sources:

- [PEP paper](https://arxiv.org/abs/2103.16745)
- [ASCL PEP listing via orbit search](https://ascl.net/code/search/orbit)
- [ILRS paper using PEP on LLR data](https://ilrs.cddis.eosdis.nasa.gov/lw19/docs/2014/Papers/3148_Martini_paper.pdf)

### Tier B: `JPL` and `POLAC/ELPN`

These are strong science-grade routes, but they do **not** look like simple
public drop-in code hand-offs in the current bounded scout.

#### `JPL`

The `JPL` Lunar Analysis Center is clearly mature:

- it fits lunar laser ranges from 1970 to the present,
- reports cm-class weighted RMS residuals,
- and distributes public ephemerides.

But the public-facing path found in this scout is:

- `DE` ephemerides,
- `SPICE`,
- and ephemeris reader tooling,

not a released end-to-end `LLR` estimator stack.

Sources:

- [JPL planetary and lunar ephemerides](https://ssd.jpl.nasa.gov/planets/eph_export.html)
- [ILRS technical report, JPL Lunar Analysis Center](https://ilrs.gsfc.nasa.gov/docs/2020/ilrsreport_2016_section7.pdf)

#### `POLAC / ELPN`

`POLAC` is also clearly mature:

- long-running `LLR` analysis center,
- `ELPN` numerical lunar ephemeris,
- variational-equation partials,
- explicit gravity-test use cases.

But in this bounded scout, no public estimator code release was identified that
could be treated as a local drop-in replacement for the frozen `MLRS` branch.

Sources:

- [ILRS technical report, POLAC / ELPN](https://ilrs.gsfc.nasa.gov/docs/2020/ilrsreport_2016_section7.pdf)
- [22nd IWLR paper using ELPN in operational prediction context](https://ilrs.gsfc.nasa.gov/lw22/papers/S10/S10-06_Bourgoin_Paper.pdf)

### Tier C: `INPOP/CALCEPH` and `EPM/IAA`

These are useful support or collaboration paths, but not the best direct
hand-off targets for the final weak-field posterior in this workspace.

#### `INPOP + CALCEPH`

What is public:

- `INPOP` ephemeris products,
- direct-download ephemeris files,
- a statement that Moon residual computation is proposed to users,
- open-source `CALCEPH` libraries for reading `INPOP` and `DE` files.

What was **not** identified:

- a released full `LLR` fit engine comparable to what `Request 4` now needs.

Sources:

- [INPOP](https://www.imcce.fr/inpop/)
- [CALCEPH intro](https://calceph.imcce.fr/docs/4.0.1/html/c/calceph.intro.html)
- [CALCEPH project page](https://www.imcce.fr/inpop/calceph)

#### `EPM / IAA RAS`

The `IAA RAS` route is mature and useful scientifically, but the official site
states plainly that:

> the scientific software engine ... is not distributed.

So this is not a local estimator-code hand-off path.

Sources:

- [IAA RAS references / online ephemeris service](https://iaaras.ru/en/dept/ephemeris/online/references/)
- [ILRS technical report, IAA RAS LLR activity](https://ilrs.gsfc.nasa.gov/docs/2020/ilrsreport_2016_section7.pdf)

## Project Consequence

The practical `Request 4` conclusion is now:

1. keep `MLRS` frozen as a bounded weak-field observation-operator laboratory,
2. use `PEP` as the first external estimator acquisition/build target,
3. if `PEP` cannot be acquired or built cleanly, switch to a collaboration /
   external-analysis-center route (`JPL` or `POLAC`),
4. treat `INPOP/CALCEPH` and `JPL DE/SPICE` as support infrastructure, not as a
   substitute for a full weak-field estimator.

## Current Status Line

The most accurate one-line status is now:

> `Request 4` should no longer grow `MLRS`; the best public-code hand-off
> target found is `PEP`, while `JPL` and `POLAC/ELPN` remain strong but
> collaboration-style external routes, and `INPOP/CALCEPH` plus `DE/SPICE`
> remain support layers rather than full posterior engines.
