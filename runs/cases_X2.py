import itertools
from _cases_X import mass_fracs, atmos

doserates = [32e5, 32e4]  # doserates
int_combos = list(itertools.product(*map(range, map(len, [mass_fracs, atmos, doserates]))))
