import itertools

scals = [
    {"loop_scaling_inter": -1.0, "loop_scaling_intra": -0.0},
]
mass_fracs = [5e-4, 20e-4, 50e-4, 200e-4] # [5e-4, 20e-4, 35e-4, 50e-4]
atmos = [
      {'O2': 0.3e-3, 'N2O': 0.0},
      {'O2': 0.0, 'N2O': 29e-3},
]

crn_combos = list(itertools.product(*map(range, map(len, [scals]))))
int_combos = list(itertools.product(*map(range, map(len, [mass_fracs, atmos]))))
