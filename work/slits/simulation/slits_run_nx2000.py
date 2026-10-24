import numpy as np
import json, os
from paraxial_reference import MODEL_VERSION, source_hash
from slits_density import build_initial_field, propagate_free, exact_reference, velocity_fields, rel_l2, fringe_spacing, run

if __name__ == '__main__':
    tag = os.environ.get('SLITS_TAG', 'nx2000')
    params = {
        'Ny': 1024, 'Ymax': 30.0, 'k': 2.0 * np.pi, 'd': 4.0, 'a': 1.0,
        'L': 20.0, 'Nx': 2000, 'Np': 100000, 'eps': 0.0, 'x_start': 0.0,
        'seed': 42,
        'check_times': [0.25, 0.5, 0.75],
    }
    from paraxial_reference import output_paths
    json_path, npz_path = output_paths(tag, __file__)
    metrics, diagnostics, issues, (y, I_num, I_ref), p = run(params)
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
    np.savez_compressed(npz_path, y=y, I_num=I_num, I_ref=I_ref)
    import hashlib
    payload['arrays'] = os.path.basename(npz_path)
    payload['arraysSha256'] = hashlib.sha256(open(npz_path, 'rb').read()).hexdigest()
    with open(json_path, 'w', encoding='utf-8') as fp:
        json.dump(payload, fp, indent=2)
    m = metrics
    print('Nx = %d, dx = %.5f' % (p['Nx'], p['dx']))
    print('prop_rel_l2        = %.3e' % m['prop_rel_l2'])
    for k2, v in m['density_errors'].items():
        print('%-20s = %s' % (k2, '%.4f' % v if v is not None else 'None'))
    print('max_density_error  = %.4f' % m['max_density_error'])
    print('n_clipped          = %d' % diagnostics['n_clipped'])
    print('max_velocity       = %.3f' % diagnostics['max_velocity'])
    print('issues             = %s' % issues)
    print('saved              = %s' % json_path)
