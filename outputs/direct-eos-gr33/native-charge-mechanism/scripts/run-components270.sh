# Phase270: component/part decomposition of the 4x endpoint charge (phase-268 T readout), same launcher as the depth tool.
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cd /home/lpaiu/work/native-refined268-runtime || exit 1
tr -d '\r' < $S/phase270-components.py > .phase270-components.py
K="baryon_g gas_nonrest_energy_erg nonrest_trace_erg nonrest_stress_erg pressure_volume_erg photon_energy_erg photon_radial_pressure_erg metric_stress_erg"
specs="all state geometry"
for k in $K inner_cumulative_energy_erg outer_cumulative_energy_erg; do specs="$specs state:$k"; done
for k in $K; do specs="$specs geometry:$k"; done
taskset -c 7 python3 .phase251-reader-launch.py .phase270-components.py readout268-quad64-work T $specs > .phase270-components.stdout.log 2> .phase270-components.stderr.log
echo "exit=$?" > .phase270-components.done
