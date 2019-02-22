#!/usr/bin/env python3

import json
import argh

from aq_radiolysis.pulsed import (
    get_t_low, get_durations_doserates_high_low
)
from chempy.units import to_unitless, patched_numpy as pnp, default_units as u

def main(
        outfile='varied.json',
        dr_pulse_Gy_s=32e6,                 # Polish samples
        t_pulse_s=5e-6,                     # Polish samples
        f_pulse_Hz=400.0,                   # Polish samples
        count_wait="[(2,0.090),(125,60)]",  # Polish samples
        tot_dose_Gy=80e3
):
    count_wait = [(int(c), float(t)*u.s) for c, t in eval(count_wait or "[]")]
    dur, dr = get_durations_doserates_high_low(
        t_pulse_s*u.s, get_t_low(t_pulse_s*u.s, f_pulse_Hz*u.Hz), dr_pulse_Gy_s*u.Gy/u.s,
        total_dose=tot_dose_Gy*u.Gy, count_wait=count_wait)
    json.dump([
        to_unitless(dur, u.s).tolist(),
        dict(doserate=to_unitless(dr, u.Gy/u.s).tolist())
    ], open(outfile, 'wt'))


if __name__ == '__main__':
    argh.dispatch_command(main)
