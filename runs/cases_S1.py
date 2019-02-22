import itertools

scals = [
    {"loop_scaling_inter": -2, "loop_scaling_intra": -2},
]
mass_fracs = [500e-6, 2000e-6, 3500e-6, 5000e-6]
atmos = [
      {'O2': 0.3e-3, 'N2O': 0.0},
      {'O2': 0.0, 'N2O': 29e-3},
]

crn_combos = list(itertools.product(*map(range, map(len, [scals]))))
int_combos = list(itertools.product(*map(range, map(len, [mass_fracs, atmos]))))
