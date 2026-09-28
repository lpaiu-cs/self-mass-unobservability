"""Phase 288: add the layer-by-layer thermal-relaxation estimate of phase 287 to Section 4.6 (draft and manuscript alike).

Four replacements, each asserted to occur exactly once in both files; the manuscript keeps its CRLF line endings.
Numbers are read from outputs/direct-eos-gr33/native-thermal-relaxation/phase287/relax.json and checked against the text.
"""
import json
from pathlib import Path
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
r = json.loads((root/'outputs/direct-eos-gr33/native-thermal-relaxation/phase287/relax.json').read_text(encoding='utf-8'))
qi, qo = r['n_in']['debye_quadrature_rel'], r['n_out']['debye_quadrature_rel']
di, do = r['n_in']['depth_km_at_tau_1_over_omega'], r['n_out']['depth_km_at_tau_1_over_omega']
mi, mo = r['n_in']['mass_fraction_above'], r['n_out']['mass_fraction_above']
assert f'{qi:.1e}' == '4.0e-09' and f'{qo:.1e}' == '3.3e-07', (qi, qo)
assert round(di, -2) == 2600 and round(do, -2) == 6600 and f'{mi:.0e}' == '2e-09' and f'{mo:.1e}' == '1.4e-07', (di, do, mi, mo)

R = [
 ("- **Conjectural.** *Deeper layers.* Thermal times grow inward to the stellar cooling time, so by continuity some deeper layer has a thermal time comparable to the orbital period. Its depth, and the strength of its coupling to \\(a_i\\), were not computed.",
  "- **Imported from prior work.** *Deeper layers.* We estimated the thermal response layer by layer. Each layer relaxes on the thermal time of the layers above it, its fully relaxed limit is bracketed by the isothermal exponent \\(P_{\\rm gas}/P\\), and the static structural response is solved again with the faster layers relaxed. The layers whose thermal time equals \\(1/\\omega\\) lie about 2,600 km deep for the inner orbit and 6,600 km for the outer one, above about \\(2\\times10^{-9}\\) and \\(1.4\\times10^{-7}\\) of the mass. Summing the shell strengths, of either sign, with the Debye weight \\(\\omega\\tau/(1+\\omega^2\\tau^2)\\) gives a lag of \\(4.0\\times10^{-9}\\,|\\mathcal S_{\\rm struct}|\\) at the inner orbital frequency and \\(3.3\\times10^{-7}\\,|\\mathcal S_{\\rm struct}|\\) at the outer one. **Conjectural.** This estimate does not use the two assumptions below. It is not a non-adiabatic calculation, and its relaxed limit can differ from thermal equilibrium by factors of order unity."),
 ("Its lagged part is bounded only under the two stated assumptions on the thermal relaxation. Deeper layers with local thermal times comparable to the orbital period are expected by continuity, but the pole structure and charge coupling of their thermal response were not computed.",
  "Its lagged part is bounded under the two stated assumptions on the thermal relaxation, and a layer-by-layer estimate puts it at \\(4.0\\times10^{-9}\\) and \\(3.3\\times10^{-7}\\) of that term at the inner and outer orbital frequencies; a non-adiabatic calculation of the deep thermal response was not made."),
 ("- a thermal relaxation strength bounded by \\(\\mathcal S_{\\rm struct}\\), with a single relaxation time;",
  "- for the lag, either a thermal relaxation strength bounded by \\(\\mathcal S_{\\rm struct}\\) with a single relaxation time, or the layer-by-layer relaxation estimate;"),
 ("notes/REQUEST244_*} through \\path{notes/REQUEST286_*} and in manifests",
  "notes/REQUEST244_*} through \\path{notes/REQUEST288_*} and in manifests"),
]
ADD_AFTER = "- photosphere: \\path{native-atmosphere-reconstruction-manifest.json}"
ADD = "- layer-by-layer thermal relaxation: \\path{native-thermal-relaxation-manifest.json}"

for path, eol in ((root/'docs/white-dwarf-free-fall-charge-section.md', '\n'), (root/'paper/manuscript.md', '\r\n')):
    t = path.read_bytes().decode('utf-8')
    for old, new in R:
        assert t.count(old) == 1, (path.name, old[:50]); t = t.replace(old, new)
    anchor = ADD_AFTER + eol
    assert t.count(anchor) == 1 and ADD not in t, path.name
    t = t.replace(anchor, anchor + ADD + eol)
    path.write_bytes(t.encode('utf-8'))
    print(path.name, 'ok')
t = (root/'docs/white-dwarf-free-fall-charge-section.md').read_text(encoding='utf-8')
t = t.replace('# Draft section (revision 6):', '# Draft section (revision 7):', 1)
t = t.replace('Status (2026-09-28): final wording after the confirmation review of revision 5 (fable5.1: accept; gpt-6-astra and opus5.5: accept after minor revision).',
              'Status (2026-09-28): revision 6 passed the confirmation review of revision 5 (fable5.1: accept; gpt-6-astra and opus5.5: accept after minor revision); revision 7 adds the layer-by-layer thermal-relaxation estimate of phase 287 after a single self-review (user decision).', 1)
t = t.replace('through `notes/REQUEST286_MANUSCRIPT_REINTEGRATION_KO.md`.', 'through `notes/REQUEST288_RELAXATION_ESTIMATE_SECTION46_KO.md`.', 1)
(root/'docs/white-dwarf-free-fall-charge-section.md').write_text(t, encoding='utf-8', newline='\n')
print('header', '(revision 7)' in t and 'REQUEST288_RELAXATION' in t)
