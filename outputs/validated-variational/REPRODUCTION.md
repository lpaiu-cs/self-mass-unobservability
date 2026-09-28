# Validated variational lever

Status: Proven. The certificate is conditional on the frozen internal GR IVP in `ivp.hex`. Each hexadecimal token is an exact serialization of a binary long double. Returned JSON interval endpoints denote binary64 values. These are initial-state derivatives at recorded epochs, not derivatives of all timing parameters or the whole timing residual.

The WSL runtime is `/home/lpaiu/work/nutimo_pilot`. CAPD source is built at `request14_capd` and the client is `request14_variational`. Build CAPD from the recorded commit:

```sh
git clone https://github.com/CAPDGroup/CAPD request14_capd
git -C request14_capd checkout 731079217a9254ea2948d742df2b170895effe7f
cmake -S request14_capd -B request14_capd/build -DCMAKE_BUILD_TYPE=Release
cmake --build request14_capd/build -j 4
```

From a fresh checkout with the same isolated runtime layout:

```sh
python3 verification/validated_variational.py build
python3 verification/validated_variational.py check
/home/lpaiu/work/nutimo_pilot/request14_variational \
  outputs/validated-variational/ivp.hex \
  outputs/validated-variational/native-rhs.hex \
  0.07 /tmp/request14-local.jsonl jacobi
python3 verification/variational_audit.py
```

`capd-build.json` contains the exact compiler arguments. To replay a historical representation, compile the source snapshot specified by `run-index.json` using its recorded build arguments, with a separate binary/output path. The index binds each run's mode, requested horizon, producer and binary. Full-span attempts request 124.162279 internal time units, slightly beyond the actual last interpolation epoch; negative tests cover the entire short backward span. Each run carries its uncertainty continuously. No restart at a midpoint is used.

Status: Imported from prior work. Runs stop with a failed full-span result if the maximum physical-coordinate Jacobian interval width exceeds 1e6, if CAPD cannot validate a remainder, or if the 1800-second run budget is reached. The broad-width threshold is a declared operational ceiling, not a scientific accuracy requirement and not a chaos diagnostic. Last successful endpoint bounds remain conditional certificates even when a later requested horizon fails.

The initial live export can be regenerated with `python3 verification/validated_variational.py export` only in a fresh isolated runtime target: it refuses to overwrite an existing `nutimo_request14_export`. The read-only export patch, native source, input hashes and build/run logs are retained. This needs the Request 13 stable engine and its existing Tempo2 installation; replaying the interval calculation from the frozen IVP does not.

Status: Proven. The `mass` mode augments the 24 state coordinates with four constant fractional-mass variables. Its local 28-dimensional Jacobian includes independent mass derivatives; it is not the Jacobian of the 28 fitted timing parameters. It uses the same baseline IVP and validated solver.

Status: Conjectural. Full D2 remains open: useful whole-span bounds, physical parameter initialization including the mapping to masses, the full delay/inverse-time map, and propagation into the nuisance projector and likelihood are still required. The manuscript PDF and submission ZIP are unchanged Request 12 artifacts.
