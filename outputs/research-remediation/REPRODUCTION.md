# Request 13 evidence

Producer environment: WSL Ubuntu 22.04, Python 3.10.12, NumPy 2.2.6. The isolated request13_deps directory contains SciPy 1.15.3 and lalsuite distribution 7.26.15 (lal 7.7.1, lalsimulation 6.2.1). The historical Python/engine environment was not replaced.

Run the lightweight artifact checks from the repository root:

    rtk proxy wsl -d Ubuntu-22.04 -- bash -lc "cd /mnt/e/lab/self-mass-unobservability && OPENBLAS_NUM_THREADS=2 python3 verification/remediation_audit.py"
    rtk proxy wsl -d Ubuntu-22.04 -- bash -lc "cd /mnt/e/lab/self-mass-unobservability && OPENBLAS_NUM_THREADS=2 python3 verification/nonlinear_joint_fit.py --check"
    rtk proxy wsl -d Ubuntu-22.04 -- bash -lc "cd /mnt/e/lab/self-mass-unobservability && OPENBLAS_NUM_THREADS=2 python3 verification/verify_unified_paper.py"

The stellar producer additionally needs PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps. Live timing producers need LD_LIBRARY_PATH to the selected isolated source directory, TEMPO2=/home/lpaiu/work/nutimo_pilot/install/third_party/tempo2 and OMP_NUM_THREADS=2. Build commands and source/library hashes are in build.json and stable-build.json. Source copies and runtime working directories are external; the patches, producer programs and residual outputs are preserved here. Builders intentionally refuse to replace existing builds.

`request12-revision-manifest.json` preserves the prior manifest. The previous manuscript PDF and source ZIP remain frozen Request 12 artifacts; they are not a completed Request 13 submission. The active paper manifest identifies the five updated supporting notes explicitly.

The default covariance range was expanded after the initial pilot reached its edge; exact zero red amplitude was then admitted. The zero_refined run used all 28 freshly recomputed derivative columns. It ended after three failed local searches, not a stationarity certificate. All producers completed; no background computation is required to obtain the stored result.

The final manifest binds raw arrays, patches, source snapshots, reports and producer code. Re-running producers changes provenance/elapsed-time logs and requires a new manifest. A successful artifact check does not turn the explicitly false derivative and physical-inference gates into successes.
