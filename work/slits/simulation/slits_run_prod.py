"""Larger corrected pilot. Uses the same audited solver as the small pilot."""
import numpy as np
import json, os
from paraxial_reference import MODEL_VERSION, source_hash
from slits_simulation import run

if __name__ == '__main__':
    tag = os.environ.get('SLITS_TAG', 'prod')
    params = {
        'Ny': 1024, 'Ymax': 30.0, 'k': 2.0 * np.pi, 'd': 4.0, 'a': 1.0,
        'L': 20.0, 'Nx': 500, 'Np': 1000000, 'eps': 0.0, 'x_start': 0.0,
        'seed': 42,
    }
    from paraxial_reference import output_paths
    json_path, npz_path = output_paths(tag, __file__)
    metrics, diagnostics, issues, (y, I_num, I_ref, H), p = run(params)
    here = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(here, 'corrected_v2_result_%s.json' % tag)
    npz_path = os.path.join(here, 'corrected_v2_arrays_%s.npz' % tag)
    payload = {
        'schema': MODEL_VERSION,
        'sourceSha256': source_hash(),
        'observable': 'paraxial fixed-x transverse density; not physical time or first arrival',
        'reference': 'Closed Gaussian superposition on infinite line',
        'params': p,
        'metrics': metrics,
        'diagnostics': diagnostics,
        'issues': issues,
    }
    if os.path.exists(json_path):
        raise FileExistsError('Use a new SLITS_TAG; preserve measured results.')
    with open(json_path, 'w') as fp:
        json.dump(payload, fp, indent=2)
    np.savez_compressed(npz_path, y=y, I_num=I_num, I_ref=I_ref, H=H)
    import hashlib
    payload['arrays'] = os.path.basename(npz_path)
    payload['arraysSha256'] = hashlib.sha256(open(npz_path, 'rb').read()).hexdigest()
    with open(json_path, 'w', encoding='utf-8') as fp:
        json.dump(payload, fp, indent=2)
    m = metrics
    print('prop_rel_l2        = %.3e' % m['prop_rel_l2'])
    print('born_rel_l2        = %.4f' % m['born_rel_l2'])
    print('fringe num         = %.4f' % m['fringe_spacing_num'])
    print('fringe ref         = %.4f' % m['fringe_spacing_ref'])
    print('fringe analytic    = %.4f' % m['fringe_spacing_analytic'])
    print('n_clipped          = %d' % diagnostics['n_clipped'])
    print('max_velocity       = %.3f' % diagnostics['max_velocity'])
    print('fraunhofer_valid   = %s' % diagnostics['fraunhofer_valid'])
    print('issues             = %s' % issues)
    print('saved              = %s' % json_path)
