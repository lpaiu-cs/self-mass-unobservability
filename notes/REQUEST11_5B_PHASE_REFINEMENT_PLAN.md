# Request 11.5b — declared refinement after the coarse phase result

The completed Request 11.5 coarse grids are preserved in phase-state-audit.json. At long lags the 24^3 three-phase grid misses maxima already sampled by the historical one-origin grid. The grids are not nested, so their maxima must not be interpreted as a monotonic domain comparison.

Register this follow-through before execution: same six lags, three fixed covariance amplitudes and unit-drive K=1 construction. For each cell use the maximizing phase from each of the three completed grids as three seeds. Perform a periodic 3x3x3 neighborhood search with initial angular step pi/12 and 17 successive halvings. At each scale allow at most 24 improving moves; retain and report any cap hit rather than silently expanding it. Include all seed values in the final maximum.

Report local-refined maxima and last-scale improvement, while retaining the analytic continuous upper bound from Request 11.5. A successful local refinement is not proof of the global maximum. No covariance/amplitude/lag retuning or new significance claim. The existing data's phase fit is descriptive; the assumed-noise interval interpretation remains conditional.
