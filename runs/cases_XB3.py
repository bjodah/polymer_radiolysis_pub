from _cases_XB import mass_fracs, atmos, prim_dat, crn_groups, crn_indices
import itertools

tpulses_s = [5e-5, 5e-4]
int_groups = [mass_fracs, atmos, tpulses_s]
int_indices = list(itertools.product(*map(range, map(len, int_groups))))
