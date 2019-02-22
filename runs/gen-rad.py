#!/usr/bin/env python3

import json
from chempy import ReactionSystem, Reaction
from chempy.kinetics.rates import MassAction, RadiolyticBase
from chempy.units import SI_base_registry, to_unitless, get_derived_unit, unitless_in_registry, default_units as u
from aq_radiolysis.models import elliot_reaction_2009 as elliot
from aq_radiolysis.util import radiolytic_yields_to_reactions
from pykinetgine.util import reac, substance, mk_kilomol_registry

def mk_rsys(*, R14a=False, R21=False, R22a=False):
    T=298.15*u.K
    rsys_elliot = ReactionSystem(*elliot.get_reactions_and_species(T=T, R14a=R14a, R21=R21, R22a=R22a))
    g_values = {k: v(T).rescale(u.mole/u.joule) for k, v in elliot.low_LET_g.items()}
    rsys_N2O = ReactionSystem.from_string("N2O + e-(aq) -> N2 + O-; 9.6e9/M/s")
    rsys = rsys_elliot + rsys_N2O + radiolytic_yields_to_reactions(g_values)
    rsys.sort_substances_inplace()
    return rsys

def dump_reactions(R14a=False, R21=False, R22a=False, output=''):
    rsys = mk_rsys(R14a=R14a, R21=R21, R22a=R22a)
    u_reg = mk_kilomol_registry()
    variables = dict(density=998*u.kg/u.m3)
    rendered = json.dumps([reac(r, u_reg, variables=variables, json=True) for r in rsys.rxns], indent=4)
    if output:
        open(output, 'wt').write(rendered)
    else:
        print(rendered)

def dump_substances(output=''):
    rsys = mk_rsys()
    rendered = json.dumps([substance(s, json=True) for s in rsys.substances.values()], indent=4)
    if output:
        open(output, 'wt').write(rendered)
    else:
        print(rendered)

if __name__ == '__main__':
    import argh
    argh.dispatch_commands([dump_reactions, dump_substances])
