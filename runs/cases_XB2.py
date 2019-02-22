from _cases_XB import mass_fracs, atmos, prim_dat, crn_groups, crn_indices
import itertools

doserates_Gy_s = [32e5, 32e4]  # doserates
int_groups = [mass_fracs, atmos, doserates_Gy_s]
int_indices = list(itertools.product(*map(range, map(len, int_groups))))
