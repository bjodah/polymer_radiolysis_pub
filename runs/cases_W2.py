import itertools

scals = [
    {"loop_scaling_inter": -1, "loop_scaling_intra": -2}, # ref
    {"loop_scaling_inter": -2, "loop_scaling_intra": -2},
    {"loop_scaling_inter": -4, "loop_scaling_intra": -2},
    {"loop_scaling_inter": -1, "loop_scaling_intra": -3},
    {"loop_scaling_inter": -1, "loop_scaling_intra": -4},
]
mass_fracs = [5e-4, 5e-3]
atmos = [
      {'O2': 0.3e-3, 'N2O': 0.0},
      {'O2': 0.0, 'N2O': 29e-3},
]

crn_combos = list(itertools.product(*map(range, map(len, [scals]))))
int_combos = list(itertools.product(*map(range, map(len, [mass_fracs, atmos]))))
