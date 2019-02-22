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
    dur, varied = json.load(open('varied-%s.json' % token))
    dr = varied.pop('doserate')
    assert not varied
    cx, cy, ip1d = mk_cumulative_interpolator(dur, dr)

    cases = __import__('cases_%s' % token)
    nscal, nmf, natmos = len(cases.scals), len(cases.mass_fracs), len(cases.atmos)
    res = lambda **kw: ('%s-res/%s_{scal}_{mf}_{atmo}.txt' % (token, token)).format(**kw)

    rs = dict()
    for iscal in range(nscal):
        for imf in range(nmf):
            for iatmo in range(natmos):
                rs[iscal, imf, iatmo] = read(res(scal=iscal, mf=imf, atmo=iatmo))

    def plot_example(iscal=0):
        fig, ax = plt.subplots(1, 1, figsize=(16, 6))
        for i, k in enumerate(rs[iscal, 0, 0]['spg']['ord_rad']):
            ax.loglog(rs[iscal, 0, 0]['tout'], rs[iscal, 0, 0]['Cout'][k], ls=['-',':','-.'][i % 3], label=k)
        ax.set_xlim([1e-10, None])
        ax.set_ylim([1e-15, 100])
        ax.legend()
        print(rs[iscal, len(cases.mass_fracs)-1, 0]['ic0'])

    def export_txt():
        for iscal in range(nscal):
            for imf in range(nmf):
                for iatmo in range(natmos):
                    r = rs[iscal, imf, iatmo]
                    tot = total_of(r['Cout'], r['tokg']['a'])
                    np.savetxt('%s_%d_pmf%d_%.1gmM-O2_%.1gmM-N2O.txt' % (
                        token,
                        iscal,
                        int(-np.log10(cases.mass_fracs[imf])),
                        cases.atmos[iatmo]['O2'],
                        cases.atmos[iatmo]['N2O'],
                    ), np.array([r['tout'], tot]).T)

    def plot_single_token(*, token, iscal=0, xlim=None, frac_remaining=None, per_number=False,
                          subplots_kw={'sharey': 'col'}):
        fig, axes = plt.subplots(nmf, natmos, figsize=(14, 14), **subplots_kw)
        if per_number and frac_remaining is not None:
            raise ValueError("Incompatibel options")
        for imf in range(nmf):
            for iatmo in range(natmos):
                ax, r = axes[imf, iatmo], rs[iscal, imf, iatmo]
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
                if iatmo == natmos-1:
                    ax.set_ylabel('mass frac = %s (%.3g M)' % (cases.mass_fracs[imf], r['ic0']))
                    ax.yaxis.set_label_coords(1.1, 0.5)
                if imf == 0:
                    ax.set_title(', '.join(['[%s] = %.1g M' % (k, cases.atmos[iatmo][k])
                                            for k in 'O2 N2O'.split()]))
                if imf == 2:
                    ax.set_xlabel('time / s')
                if iatmo == 0:
                    ax.set_ylabel('concentration / M')
                if xlim:
                    ax.set_xlim(xlim)

    def plot_tokens(iscal=0):
        fig, all_axes = plt.subplots(len(tokens)-1, natmos, figsize=(16, 3*len(tokens)))
        for tki, (axes, tk) in enumerate(zip(all_axes, tokens[1:])):
            for iatmo, ax in zip(range(natmos), axes):
                for imf in range(nmf):
                    r = rs[iscal, imf, iatmo]
                    c = viridis(imf/(nmf-1.0))
                    ax.plot(r['tout'], total_of(r['Cout'], r['tokg'][tk])/r['ic0'],
                            label='%.1g ppm polymer' % (cases.mass_fracs[imf] * 1e6), c=c)
                ax.hlines(sorted(map(float, r['tokg'][tk].keys())), 0, 1,
                          transform=ax.get_yaxis_transform(),
                          linestyle='--', alpha=0.5, color='k', lw=0.5)
                ax.legend()
                if tki == 0:
                    ax.set_title(', '.join(['[%s] = %.1g M' % (k, cases.atmos[iatmo][k])
                                            for k in 'O2 N2O'.split()]))
                if tki == len(tokens)-2:
                    ax.set_xlabel('time / s')
                if iatmo == 0:
                    ax.set_ylabel(tk)

                ax.set_yscale('log')
                ax.set_ylim([1e-12, 1e2])

    def plot_conc(sks={'H2O2': (1e-6, 1e-3, 'log'), 'O2': (1e-15, 1e-4, 'log'),
                       'N2O': (28.95e-3, 29.01e-3, 'linear'),
                       'HO2': (1e-12, 1e-5, 'log'), 'OH': (1e-16, 1e-4, 'log'), 'H+': (1e-10, 1e-4, 'log')},
                      iscal=0):
        fig, all_axes = plt.subplots(len(sks)+1, natmos, figsize=(16, 3*len(sks)))
        for iatmo in range(natmos):
            for imf in range(nmf):
                r = rs[iscal, imf, iatmo]
                lbl = 'mass frac: %s' % cases.mass_fracs[imf]
                c = viridis(imf/(nmf-1.0))
                for ski, (sk, (lb, ub, yscl)) in enumerate(sks.items()):
                    ax = all_axes[ski, iatmo]
                    ax.plot(r['tout'], r['Cout'][sk], label=lbl, c=c)
                    ax.set_ylim([lb, ub])
                    ax.set_yscale(yscl)
                    ax.set_xlim([1e-9, 300])
                    ax.set_xscale('log')
                    if iatmo == 0 and imf == 0:
                        ax.set_ylabel('[%s] / M' % sk)
                all_axes[-1, iatmo].plot(r['tout'], average_length_by_weight(r['Cout'], r['tokg']),
                                        label=lbl, c=c)
            all_axes[-1, iatmo].set_xlabel('time / s')
            all_axes[0, iatmo].set_title(', '.join(['[%s] = %.1g M' % (k, cases.atmos[iatmo][k])
                                            for k in 'O2 N2O'.split()]))
            all_axes[-1, iatmo].axhline(int(r['sk0'][1:].split('a')[0]), ls='--', lw=0.5,
                                        alpha=0.5, c='k')
            #all_axes[-1, iatmo].set_ylim([0, 500])
            all_axes[-1, iatmo].set_yscale('log')
        all_axes[-1, 0].set_ylabel('Average length by weight')
        _ = all_axes[-1, 0].legend()

    # This cell generates figures for Mats' conference presentation
    def plot_Mats(iscal=0, zoom=False, xlim_zoom=[0, 10e-3], xlim_kGy=None,
                  species={'H2O2': 'H_2O_2', 'e-(aq)': 'e^-(aq)', 'N2O': 'N_2O', 'O2': 'O_2', 'H+': 'H^+', 'OH': '^.OH'},
                  ylim={'O2': [1e-10, 1e-4], 'OH': [0, 2e-5]}, yscale={'OH': 'linear'}):
        for imf in range(nmf):
            if zoom:
                if imf != 0:
                    # N2O, lägsta mf, 90ms - 110 ms
                    continue
                ext = '_zoom.png'
            else:
                ext = '.png'

            ppm = cases.mass_fracs[imf] * 1e6
            if cases.atmos[0]['O2'] == 0:
                initro, iair = 0, 1
                assert cases.atmos[iair]['N2O'] == 0
            else:
                iair, initro = 0, 1
                assert cases.atmos[initro]['O2'] == 0
            air, nit = rs[iscal, imf, iair], rs[iscal, imf, initro]
            fig1, ax1 = plt.subplots(1, 1)
            for lbl, cont in list(zip(['air', 'nitrous oxide'], [air, nit])):
                n0 = int(cont['sk0'][1:].split('a')[0])  # initial length
                ax1.axhline(n0, ls='--', lw=0.5, alpha=1.0 if zoom else 0.5, c='k')
                ax1.plot(cont['tout'] if zoom else ip1d(cont['tout'])/1e3, average_length_by_weight(cont['Cout'], cont['tokg']), label=lbl)
            ax1.set_ylabel('Chain length / number of segments')
            ax1.legend()
            if zoom:
                ax1.set_xlabel('Time / s')
                ax1.set_xlim(xlim_zoom)
            else:
                ax1.set_xlabel('Dose / kGy')
                ax1.set_xlim(xlim_kGy)
            fig1.savefig('poly_%s_size_%dppm%s' % (token, ppm, ext), dpi=300)

            for specie, sp_tex in species.items():  # ['H2O2', 'e-(aq)', 'N2O', 'O2', 'H+']:
                fig2, ax2 = plt.subplots(1, 1)
                for lbl, cont in list(zip(['air', 'nitrous oxide'], [air, nit])):
                    ax2.plot(cont['tout'] if zoom else ip1d(cont['tout'])/1e3, cont['Cout'][specie], label=lbl, alpha=1.0 if zoom else 0.7)
                ax2.set_ylabel(r'[$\mathrm{%s}$] / M' % sp_tex)
                ax2.ticklabel_format(axis='y', style='sci', scilimits=(-4, 4))
                ax2.set_yscale((yscale or {}).get(specie, 'log'))
                ax2.set_ylim(ylim.get(specie, [1e-8, 3e-2]))
                ax2.legend()
                fig2.tight_layout()
                if zoom:
                    ax2.set_xlim(xlim_zoom)
                    ax2.set_xlabel('Time / s')
                else:
                    ax2.set_xlim(xlim_kGy)
                    ax2.set_xlabel('Dose / kGy')

                fig2.savefig('poly_%s_%s_%dppm%s' % (token, specie, ppm, ext), dpi=300)

            fig3, ax3 = plt.subplots(1, 1)
            for lbl, cont in list(zip(['air', 'nitrous oxide'], [air, nit])):
                x, y = cont['tout'] if zoom else ip1d(cont['tout'])/1e3, sum_up(cont['Cout'], cont['ctrg']['unsaturated'])/cont['ic0']
                ax3.plot(x, y, label=lbl)
                if not zoom:
                    p = np.polyfit(x, y, 1)
                    ax3.plot(x[[0,-1]], np.polyval(p, x[[0,-1]]), c='k', label="G=%.5g mol/J" % (p[0]/1e3/0.998*cont['ic0']))
            ax3.set_ylabel('Number of double bonds per polymer')
            ax3.legend()
            if zoom:
                ax3.set_xlim(xlim_zoom)
                ax3.set_xlabel('Time / s')
            else:
                ax3.set_xlim(xlim_kGy)
                ax3.set_xlabel('Dose / kGy')
            fig3.savefig('poly_%s_unsaturated_%dppm%s' % (token, ppm, ext), dpi=300)

            fig4, ax4 = plt.subplots(1, 1)
            if zoom:
                # dose vs time
                ax4.plot(cont['tout'], ip1d(cont['tout'])/1e3)
                ax4.set_ylabel('Dose / kGy')
                ax4.set_xlabel('time / s')
                ax4.set_xlim(xlim_zoom)
                fig4.savefig('dose%s' % ext, dpi=300)
            else:
                # final size distribution:
                for lbl, cont in list(zip(['air', 'nitrous oxide'], [air, nit])):
                    nrm = 0.0
                    arr = []
                    for w, g in cont['tokg']['n'].items():
                        val = 0.0
                        for k in g:
                            val += w*w*cont['Cout'][k][-1]
                            nrm += w*cont['Cout'][k][-1]
                        arr.append([w, val])
                    arr = np.array(arr)
                    n0 = int(cont['sk0'][1:].split('a')[0])  # initial length
                    ax4.plot(arr[:, 0]/n0, arr[:, 1], label=lbl, marker='o')
                ax4.set_title('%d ppm' % ppm)
                ax4.legend()
                ax4.set_xlabel("Relative chain-length")
                ax4.set_ylabel("Fraction of specific length")
                fig4.savefig('poly_%s_final_size_distribution_%dppm%s' % (token, ppm, ext), dpi=300)

            for lbl, cont, colr in list(zip(['air', 'nitrous oxide'], [air, nit], 'tab:blue tab:orange'.split())):
                fig5, ax5 = plt.subplots(1, 1)

                tot = total_of(cont['Cout'], cont['tokg']['a'])
                ax5.plot(cont['tout'] if zoom else ip1d(cont['tout'])/1e3, tot/cont['ic0'], label=lbl, alpha=1.0 if zoom else .5, lw=1.0 if zoom else 0.5, c=colr)

                ax5.set_ylabel('Average number of alkyl radicals per chain')
                ax5.legend()
                if zoom:
                    ax5.set_xlim(xlim_zoom)
                    ax5.set_xlabel('Time / s')
                else:
                    ax5.set_xlim(xlim_kGy)
                    ax5.set_xlabel('Dose / kGy')
                fig5.savefig('poly_%s_alkyl_%s_%dppm%s' % (token, lbl, ppm, ext), dpi=300)

    return AttributeContainer(**dict([(k, v) for k, v in locals().items() if not k.startswith('_')]))
