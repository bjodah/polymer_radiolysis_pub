import itertools
from _cases_X import mass_fracs, atmos

tpulses = [5e-5, 5e-4]  # t_pulse_s
int_combos = list(itertools.product(*map(range, map(len, [mass_fracs, atmos, tpulses]))))
