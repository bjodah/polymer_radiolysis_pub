#!/usr/bin/env python3

import json
import matplotlib
import matplotlib.pyplot as plt
import argh
from analysis import parse_results, get_token_groups, animate, save_mov


def main(results="results.txt", savemov="None", fps=30, dpi=100, text_conc="N2O,O2,H2O2,HO2,O2-,OH,e-(aq),H",
         vmin=1e-16, vmax=1e-5, varied='varied.json', interactive=False, delta=False):
    tout, Cout, info = parse_results(results)
    tokg, ctrg, spg = get_token_groups(Cout)
    kw = dict(
        text_conc=[_ for _ in text_conc.split(',') if _ != ""],
        varied=json.load(open(varied)),
    )
    if delta:
        kw['yscale'] = 'symlog', dict(linthreshy=1e-10)
        kw['conc_bounds'] = [-1e-3, 1e-3]
        kw['cmap'] = 'coolwarm'
        kw['norm'] = matplotlib.colors.SymLogNorm(linthresh=vmin, vmin=-vmax, vmax=vmax)
    else:
        kw['yscale'] = 'log', {}
        kw['conc_bounds'] = [vmin, 30e-3]
        kw['cmap'] = 'viridis'
        kw['norm'] = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)

    if savemov == 'None':
        anim1 = animate(tout, Cout, info, tokg, interactive=interactive, **kw)
        plt.show()
    else:
        save_mov(tout, Cout, info, tokg, savemov=savemov, fps=fps, dpi=dpi, **kw)

if __name__ == '__main__':
    argh.dispatch_command(main)
