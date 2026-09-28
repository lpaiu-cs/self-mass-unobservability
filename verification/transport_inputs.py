"""Source-bound transport input contracts; not a native intervention result."""
import json, shutil, sys
import numpy as np
import direct_eos_gr as g


def audit():
    out=g.OUT/'transport-inputs';out.mkdir(parents=True,exist_ok=True)
    files=['kap/public/kap_lib.f90','kap/private/kap_eval.f90',
        'star/private/opacities.f90','neu/public/neu_lib.f90','star/private/neu.f90']
    hashes={}
    for rel in files:
        source=g.c.fresh.MESA/rel;target=out/'sources'/rel;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(source,target);hashes[rel]=g.c.sha(target)
    kap=(out/'sources/kap/public/kap_lib.f90').read_text()
    neu=(out/'sources/neu/public/neu_lib.f90').read_text()
    opacity=kap.split('subroutine kap_get_Type1(',1)[1].split('end subroutine kap_get_Type1',1)[0]
    assert 'free_e := total combined number per nucleon of free electrons and positrons' in opacity
    assert all(name in opacity for name in ['lnfree_e','d_lnfree_e_dlnRho','d_lnfree_e_dlnT','zbar'])
    signature=neu.split('subroutine neu_get(',1)[1].split(')',1)[0]
    assert all(name in signature for name in ['abar','zbar','z2bar','log10_Tlim'])
    assert 'eta' not in signature and 'free_e' not in signature
    assert 'for T < 10^7, the neutrino losses are simply set to 0' in neu
    assert 'kap = 1d0 / (1d0/kap_rad + 1d0/kap_ec)' in (out/'sources/kap/private/kap_eval.f90').read_text()
    data=dict(np.load(g.OLD/'initial-state.npz'));below=data['lnT']<np.log(1e7)
    result=dict(classification='Imported from prior work',source_contract_checked=True,
        source_sha256=hashes,state_sha256=g.c.sha(g.OLD/'initial-state.npz'),
        opacity='The source Type1/Type2 blend requires total free electron PLUS positron number per nucleon and its log derivatives; zbar also enters electron conduction. Nuclear net Ye or eta alone is not this input contract.',
        neutrino='The source neu_get takes T,rho,abar,zbar,z2bar, a temperature cutoff and type flags. It has no external eta/free-electron input. Its documented low-temperature zero is a cutoff of the declared fit, not a physical zero-loss theorem.',
        cells_below_documented_neutrino_fit_temperature=int(below.sum()),
        baryon_fraction_below_documented_neutrino_fit_temperature=float(data['dm'][below].sum()/data['dm'].sum()),
        total_opacity_includes_conduction=True,
        actual_new_EOS_opacity_inputs_intervened=False,physical_opacity_or_neutrino_errors_certified=False,
        next_boundary='Bind the actual opacity call and its new EOS electron/positron inputs before a nonlinear transport claim. Keep the thermal-neutrino fit assumptions and missing low-temperature physical error separate. Do not double-count conduction.')
    (out/'audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print('TRANSPORT INPUT CONTRACTS',result['cells_below_documented_neutrino_fit_temperature'],
        'cells below documented neutrino fit cutoff',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
