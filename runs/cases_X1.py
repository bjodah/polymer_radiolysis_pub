import itertools
from _cases_X import mass_fracs, atmos

freqs = [400, 40, 4]  # frequencies in Hz
int_combos = list(itertools.product(*map(range, map(len, [mass_fracs, atmos, freqs]))))
