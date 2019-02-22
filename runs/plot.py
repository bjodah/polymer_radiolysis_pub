#!/usr/bin/env python3

import argh
import numpy as np
import matplotlib.pyplot as plt


def main(results='results.txt', ls='-,--,:', # ,-.
         c=('tab:blue,tab:orange,tab:green,tab:red,tab:purple,'
            'tab:brown,tab:pink,tab:gray,tab:olive,tab:cyan,black'),
            savefig='None', dpi=100,
             text_conc="H,H+,H2,H2O,H2O2,HO2,HO2-,N2,N2O,O-,O2,O2-,O3-,OH,OH-,e-(aq)"):
    text_conc=text_conc.split(',')
    ls = ls.split(','); c=c.split(',')
    with open(results) as fh:
        keys = fh.readline().lstrip('#').split()[1:]  # drop t
    data = np.genfromtxt(results)
    t = data[:, 0]
    ally = dict(zip(keys, data[:, 1:].T))
    fig, ax = plt.subplots(1, 1, figsize=(16, 10))
    i = 0
    for k, y in ally.items():
        if text_conc and k not in text_conc:
            continue
        if np.any(y > 0):
            ax.plot(t, y, label=k, c=c[i%len(c)], ls=ls[i%len(ls)])
            i += 1
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    if savefig == 'None':
        fig.show()
        input('Press any key to exit')
    else:
        fig.savefig(savefig, dpi=dpi)


if __name__ == '__main__':
    argh.dispatch_command(main)
