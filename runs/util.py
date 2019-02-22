import json
from chempy import ReactionSystem, Reaction
from chempy.kinetics.rates import MassAction, RadiolyticBase
from chempy.units import SI_base_registry, to_unitless, get_derived_unit, unitless_in_registry, default_units as u
from aq_radiolysis.models import elliot_reaction_2009 as elliot
from aq_radiolysis.util import radiolytic_yields_to_reactions
from pykinetgine.util import reac, substance

def mk_rsys():
    T=298.15*u.K
    rsys_elliot = ReactionSystem(*elliot.get_reactions_and_species(T=T))
    g_values = {k: v(T).rescale(u.mole/u.joule) for k, v in elliot.low_LET_g.items()}
    rsys_elliot += radiolytic_yields_to_reactions(g_values)
    return rsys_elliot

def dump_reactions():
    rsys = mk_rsys()
    json.dump([reac(r, json=True) for r in rsys.rxns],
              open('reactions.json', 'wt'), indent=4)

def dump_substances():
    rsys = mk_rsys()
    json.dump([substance(s, json=True) for s in rsys.substances.values()],
              open('substances.json', 'wt'), indent=4)


if __name__ == '__main__':
    import argh
    argh.dispatch_commands([dump_reactions, dump_substances])
