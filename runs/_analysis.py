import json
import os
from collections import defaultdict
from IPython.display import HTML
import numpy as np
import matplotlib
from matplotlib import animation
import matplotlib.pyplot as plt
from matplotlib.cm import viridis, inferno

from chempy.util.pyutil import AttributeContainer

from analysis import (
    parse_results, get_token_groups, plot_token_groups, animate, average_length_by_weight,
    tokens, counter_tokens, total_of, sum_up, mk_cumulative_interpolator
)

def read(path):
    tout, Cout, info = parse_results(path)
    tokg, ctrg, spg = get_token_groups(Cout)
    init = []
    for sk in spg['polymeric']:
        if Cout[sk][0] > 0:
            init.append((sk, Cout[sk][0]))
    (sk0, ic0), = init
    return locals()

def make_context(token):
    cases = __import__('cases_%s' % token)
    nmf, npd = len(cases.mass_fracs), len(cases.pulse_durs)
    res = lambda **kw: ('%s-res/%s_{mf}_{tp}.txt' % (token, token)).format(**kw)

    rs = dict()
    for mf in range(nmf):
        for tp in range(npd):
            rs[mf, tp] = read(res(mf=mf, tp=tp))

    def plot_example():
        fig, ax = plt.subplots(1, 1, figsize=(16, 6))
        for i, k in enumerate(rs[0, 0]['spg']['ord_rad']):
            ax.loglog(rs[0, 0]['tout'], rs[0, 0]['Cout'][k], ls=['-',':','-.'][i % 3], label=k)
        ax.set_xlim([1e-10, None])
        ax.set_ylim([1e-15, 100])
        ax.legend()
        print(rs[len(cases.mass_fracs)-1, 0]['ic0'])

    def export_txt():
        for mf in range(nmf):
            for tp in range(npd):
                r = rs[mf, tp]
                tot = total_of(r['Cout'], r['tokg']['a'])
                np.savetxt('%s_pmf%d_pt%d.txt' % (
                    token,
                    int(-np.log10(cases.mass_fracs[mf])),
                    int(-np.log10(cases.pulse_durs[tp]))
                ), np.array([r['tout'], tot]).T)

    def plot_single_token(*, token, xlim=None, frac_remaining=None, per_number=False,
                          subplots_kw={'sharey': 'col'}):
        fig, axes = plt.subplots(nmf, npd, figsize=(14, 14), **subplots_kw)
        if per_number and frac_remaining is not None:
            raise ValueError("Incompatibel options")
        for mf in range(nmf):
            for tp in range(npd):
                ax, r = axes[mf, tp], rs[mf, tp]
                if per_number:
                    numbers = sorted(r['tokg'][token].keys())
                    for ni in range(len(numbers)):
                        _y = np.zeros(r['tout'].size)
                        for sk in r['tokg'][token][numbers[ni]]:
                            _y += r['Cout'][sk]*numbers[ni]
                        ax.plot(r['tout'], _y, label='%d * [%s: %d]' % (numbers[ni], token, numbers[ni]),
                                c=inferno(ni/(len(numbers) - 1.0)))
                    ax.set_facecolor('xkcd:cloudy blue')
                    ax.legend()
                else:
                    x, y = r['tout'], total_of(r['Cout'], r['tokg'][token])
                    until_i = -1
                    if frac_remaining is not None:
                        mi = np.argmax(y)
                        until_y = y[mi]*frac_remaining
                        until_b = y[mi:] < until_y
                        if np.any(until_b):
                            until_i = mi + np.argmax(until_b)
                    ax.plot(x[:until_i], y[:until_i])
                ax.ticklabel_format(axis='x', style='sci', scilimits=(-4, 4))
                ax.ticklabel_format(axis='y', style='sci', scilimits=(-4, 4))
                if tp == 2:
                    ax.set_ylabel('mass frac = %s (%.3g M)' % (cases.mass_fracs[mf], r['ic0']))
                    ax.yaxis.set_label_coords(1.1, 0.5)
                if mf == 0:
                    ax.set_title('t pulse = %s' % str(cases.pulse_durs[tp]))
                if mf == 2:
                    ax.set_xlabel('time / s')
                if tp == 0:
                    ax.set_ylabel('concentration / M')
                if xlim:
                    ax.set_xlim(xlim)

    def plot_tokens():
        fig, all_axes = plt.subplots(len(tokens)-1, npd, figsize=(16, 3*len(tokens)))
        for tki, (axes, tk) in enumerate(zip(all_axes, tokens[1:])):
            for pdi, ax in zip(range(npd), axes):
                for mfi in range(nmf):
                    r = rs[mfi, pdi]
                    c = viridis(mfi/(nmf-1.0))
                    ax.plot(r['tout'], total_of(r['Cout'], r['tokg'][tk])/r['ic0'],
                            label='%s ppm polymer' % (cases.mass_fracs[mfi] * 1e6), c=c)
                ax.hlines(sorted(map(float, r['tokg'][tk].keys())), 0, 1,
                          transform=ax.get_yaxis_transform(),
                          linestyle='--', alpha=0.5, color='k', lw=0.5)
                ax.legend()
                if tki == 0:
                    ax.set_title('Pulse duration: %s s' % cases.pulse_durs[pdi])
                if tki == len(tokens)-2:
                    ax.set_xlabel('time / s')
                if pdi == 0:
                    ax.set_ylabel(tk)

                ax.set_yscale('log')
                ax.set_ylim([1e-12, 1e2])

    def plot_conc(sks={'H2O2': (1e-7, 1e-4, 'log'), 'O2': (1e-17, 1e-5, 'log'),
                       'N2O': (28.95e-3, 29.01e-3, 'linear'),
                       'HO2': (1e-16, 1e-7, 'log'), 'OH': (1e-16, 1e-4, 'log')}):
        fig, all_axes = plt.subplots(len(sks)+1, npd, figsize=(16, 3*len(sks)))
        for pdi in range(npd):
            for mfi in range(nmf):
                r = rs[mfi, pdi]
                lbl = 'mass frac: %s' % cases.mass_fracs[mfi]
                c = viridis(mfi/(nmf-1.0))
                for ski, (sk, (lb, ub, yscl)) in enumerate(sks.items()):
                    ax = all_axes[ski, pdi]
                    ax.plot(r['tout'], r['Cout'][sk], label=lbl, c=c)
                    ax.set_ylim([lb, ub])
                    ax.set_yscale(yscl)
                    ax.set_xlim([1e-9, 300])
                    ax.set_xscale('log')
                    if pdi == 0 and mfi == 0:
                        ax.set_ylabel('[%s] / M' % sk)
                all_axes[-1, pdi].plot(r['tout'], average_length_by_weight(r['Cout'], r['tokg']),
                                       label=lbl, c=c)
            all_axes[-1, pdi].set_xlabel('time / s')
            all_axes[0, pdi].set_title('Pulse duration: %s s' % cases.pulse_durs[pdi])
            all_axes[-1, pdi].axhline(int(r['sk0'][1:].split('a')[0]), ls='--', lw=0.5,
                                      alpha=0.5, c='k')
            #all_axes[-1, pdi].set_ylim([0, 500])
            all_axes[-1, pdi].set_yscale('log')
        all_axes[-1, 0].set_ylabel('Average length by weight')
        _ = all_axes[-1, 0].legend()

    return AttributeContainer(**dict([(k, v) for k, v in locals().items() if not k.startswith('_')]))
