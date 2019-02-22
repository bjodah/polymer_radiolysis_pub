#!/usr/bin/env python3

import json
import argh

rad_keys = 'H H+ H2 H2O H2O2 HO2 HO2- N2 N2O O- O2 O2- O3- OH OH- e-(aq)'.split()

def render_abstol(
        outfile='abstol.json', n_levels="[1,10,100]", n0=10,
        counter_atol=1e-12, rad_atol=1e-19, base_poly_atol=1e-14,
        atol={"H2O": 55e-6, "H+": 1e-13, "OH-": 1e-13,
              "H2O2": 1e-13, "HO2-": 1e-19,
              "H2": 1e-13,
              "O2": 1e-14, "O2-": 1e-15, "HO2": 1e-15},
        default=1e-23
):
    if isinstance(n_levels, str):
        n_levels = eval(n_levels)

    if not n0 in n_levels:
        raise ValueError("Impossible")
    atol = atol.copy()
    for k in rad_keys:
        if k not in atol:
            atol[k] = rad_atol
    atol["n%da0p0o0l0" % n0] = base_poly_atol
    for n in n_levels:
        atol['tiny%d' % n] = counter_atol
        atol['scission_loop%d' % n] = counter_atol
        atol['damage_fivering%d' % n] = counter_atol
        atol['hydroxy%d' % n] = counter_atol
        atol['keto%d' % n] = counter_atol
        atol['unsaturated%d' % n] = counter_atol
    json.dump([atol, default], open(outfile, 'wt'))

if __name__ == '__main__':
    argh.dispatch_command(render_abstol)
