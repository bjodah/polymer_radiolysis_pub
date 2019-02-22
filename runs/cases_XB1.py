from _cases_XB import mass_fracs, atmos, prim_dat, crn_groups, crn_indices
import itertools

freqs = [400, 40, 4]  # frequencies in Hz
int_groups = [mass_fracs, atmos, freqs]
int_indices = list(itertools.product(*map(range, map(len, int_groups))))
