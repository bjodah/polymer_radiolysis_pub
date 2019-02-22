import os
from collections import defaultdict

import numpy as np
from scipy.interpolate import interp1d
import matplotlib
from matplotlib import animation as manimation
import matplotlib.pyplot as plt

_CI = os.environ.get('CI', '')
try:
    if len(_CI) >= 1 and _CI.lower()[1] in ('t', '1'):
        raise ImportError  # don't litter CI log with progress bar
    from tqdm import trange
except ImportError:
    trange = range


tokens = 'napol'
counter_tokens = 'tiny scission_loop damage_fivering hydroxy keto unsaturated'.split()

def parse_results(path):
    with open(path) as ifh:
        keys = ifh.readline().lstrip('#').split()[1:]
        info = ifh.readline()[1:].replace(' ', '').replace(',,', ',')
        info = eval('dict(%s)' % info)
    r = np.genfromtxt(path)
    tout = r[:, 0]
    Cout = {k: r[:, i] for i, k in enumerate(keys, start=1)}
    return tout, Cout, info

def get_token_groups(keys):
    polymeric = set()
    tokg = {t: defaultdict(list) for t in tokens}
    ctrg = {c: defaultdict(list) for c in counter_tokens}
    counters = set()
    for k in keys:
        if k.startswith('n'):
            polymeric.add(k)
            s = k[1:]
            for i in range(1, len(tokens)):
                cur, s = s.split(tokens[i])
                tokg[tokens[i-1]][int(cur)].append(k)
            tokg[tokens[-1]][int(s)].append(k)
        else:
            for beg in counter_tokens:
                if k.startswith(beg):
                    counters.add(k)
                    ctrg[beg][int(k.split(beg)[1])].append(k)
    # wg = defaultdict(list)
    # for k in polymeric:
    #     wg[int(k[1:].split('a')[0])].append(k)
    non_polymeric = set(keys) - polymeric
    ord_rad = non_polymeric - counters # ordinary radical species from water radiolysis
    return tokg, ctrg, dict(
        polymeric=polymeric,
        counters=counters,
        ord_rad=ord_rad
    )


def plot_token_groups(tout, Cout, tokg):
    fig, axes = plt.subplots(1, len(tokens), figsize=(16, 6))
    for t, ax in zip(tokens, axes):
        for k, l in tokg[t].items():
            tot = np.zeros_like(tout)
            for sk in l:
                tot += Cout[sk]
            ax.plot(tout, tot, label=str(k))
        ax.set_title(t)
        ax.legend()
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim([1e-12, 1.0])
    return fig, axes


def mk_cumulative_interpolator(xspans, dydx):
    if len(xspans) != len(dydx):
        raise ValueError("Arguments need to have equal length")
    cx = np.cumsum(xspans)
    cx = np.concatenate((cx[:1]*0, cx))
    cy = np.cumsum(np.asanyarray(xspans)*np.asanyarray(dydx))
    cy = np.concatenate((cy[:1]*0, cy))
    ip1d = interp1d(cx, cy, fill_value='extrapolate')
    return cx, cy, ip1d

c = ('tab:blue,tab:orange,tab:green,tab:red,tab:purple,'
     'tab:brown,tab:pink,tab:gray,tab:olive,tab:cyan,black').split(',')
ls = '-,--,:'.split(',')

def mpl_outside_legend(ax, **kwargs):
    """ Places a legend box outside a matplotlib Axes instance. """
    box = ax.get_position()
    ax.set_position([box.x0 + box.width*0.75, box.y0, box.width * 0.75, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='upper left', bbox_to_anchor=(-5*box.width, 1), **kwargs)

def projected_analysis(tout, Cout, info, tokg):
    R = {}
    for kix, kx in enumerate(tokens[:-1]):
        for kiy, ky in enumerate(tokens[kix+1:], kix):
            xs = sorted(tokg[kx])
            ys = sorted(tokg[ky])
            R[kx, ky] = np.zeros((tout.size, len(ys), len(xs)))
            for ix, x in enumerate(xs):
                for iy, y in enumerate(ys):
                    for k in Cout:
                        if k in tokg[kx][x] and k in tokg[ky][y]:
                            R[kx, ky][:, iy, ix] += Cout[k]
    return R


def anim_token_groups(tout, Cout, info, tokg, *, text_conc=(), varied, norm, conc_bounds,
                      unit_txt=dict(doserate='Gy/s'), cmap='viridis', yscale=('log', {})):
    fig, axes = plt.subplots(len(tokens)-1, len(tokens)-1, figsize=(14, 14))
    R = projected_analysis(tout, Cout, info, tokg)
    im = {}
    for kix, kx in enumerate(tokens[:-1]):
        for kiy, ky in enumerate(tokens[kix+1:], kix):
            xs = sorted(tokg[kx])
            ys = sorted(tokg[ky])
            nx = len(xs)
            ny = len(ys)
            im[kx, ky] = axes[kiy, kix].imshow(R[kx, ky][0, ...], extent=(0, nx, 0, ny),
                                               norm=norm, cmap=cmap)
            axes[kiy, kix].set_xlabel(kx.upper())
            axes[kiy, kix].set_ylabel(ky.upper())
            axes[kiy, kix].set_xticks(np.arange(nx)+0.5)
            axes[kiy, kix].set_xticklabels(xs)
            axes[kiy, kix].set_yticks(np.arange(ny)+0.5)
            axes[kiy, kix].set_yticklabels(ys[::-1])
            axes[kiy, kix].set_facecolor('k')
            plt.colorbar(axes[kiy, kix].get_children()[-2], ax=axes[kiy, kix])

    traces = []
    if text_conc:
        for ax in axes[:2, -1]:
            for i, k in enumerate(text_conc):
                ax.plot(tout, Cout[k], color=c[i%len(c)], ls=ls[i%len(ls)], lw=1, label=k)
            traces.extend(ax.plot([0, 0], conc_bounds, c='k', ls='--'))
            ax.set_xlim([tout[1], tout[-1]])
            ax.set_ylim(conc_bounds)
            ax.set_xlabel('Time / s')
            ax.set_ylabel('Concentration / M')
            ax.set_yscale(yscale[0], **yscale[1])
        axes[0, -1].set_xscale('log')
        mpl_outside_legend(axes[0, -1])

    for kix in range(1, len(tokens)-1):
        for kiy in range(kix):
            axes[kiy, kix].set_visible(kix == len(tokens) - 2 and kiy <= 2)
    ttl_tmplate = 't = {:10.4e}, D = {:10.4e}'

    durations, varied_values = varied
    durations = np.array(durations)
    varied_values = {k: np.array(v) for k, v in varied_values.items()}
    cx, cy, ip1d = mk_cumulative_interpolator(durations, varied_values['doserate'])
    Dtot = ip1d(tout)
    if info['success']:
        if not Dtot.size == info['npoints']*durations.size + 1:
            raise ValueError("Output size mismatch, sure you're using the right varied?")
        if (tout.size - 1) % len(varied[0]) != 0:
            raise ValueError("Bug?")

    vlines = []
    vlines.extend(axes[2, -1].plot(cx, cy*1e-3))
    vlines.extend(axes[2, -1].plot([0.0, 0.0], [0, cy[-1]*1e-3], c='k', ls='--', lw=1))
    axes[2, -1].set_xlabel("Time / s")
    axes[2, -1].set_ylabel("Dose / kGy")

    def _var_txt(i):
        return '\n'.join('lg({:>8} / {}): {:.3f}'.format(k, unit_txt.get(k, 'a.u.'), np.log10(v[i//info['npoints']]))
                         for k, v in varied_values.items())

    def _conc_txt(i):
        return '\n'.join('lg({:>8} /   M): {:.3f}'.format(k, np.log10(np.clip(Cout[k][i], 0, None))) for k in text_conc)

    s = _var_txt(0) + '\n'
    if text_conc:
        s += _conc_txt(0)
    txt = fig.text(.3, .75, s, fontsize=15, family='monospace')
    fig.tight_layout()
    def update_surf(i):
        for kix, kx in enumerate(tokens):
            for kiy, ky in enumerate(tokens[kix+1:], kix):
                im[kx, ky].set_data(R[kx, ky][i, ...])
        s = ttl_tmplate.format(tout[i], Dtot[i]) + '\n\n'
        s += _var_txt(i-1) + '\n'
        if text_conc:
            s += _conc_txt(i)

        vlines[-1].set_data([tout[i], tout[i]], [0, cy[-1]*1e-3])
        vlines[0].set_data(cx, cy*1e-3)
        for trace in traces:
            trace.set_data([tout[i], tout[i]], conc_bounds)

        txt.set_text(s)
        return list(im.values()) + traces + vlines

    return fig, update_surf


def animate(tout, Cout, info, tokg, interactive=False, **kwargs):
    kw = dict(
        frames=tout.size, interval=1, blit=False
    )
    if interactive:
        from _player import Player as Cls
        kw['maxi'] = tout.size - 1
    else:
        Cls = manimation.FuncAnimation
    return Cls(*anim_token_groups(tout, Cout, info, tokg, **kwargs), **kw)


def save_mov(tout, Cout, info, tokg, *, savemov, fps, dpi, **kwargs):
    fig, update = anim_token_groups(tout, Cout, info, tokg, **kwargs)
    metadata = dict(title='%s - %s' % (savemov, os.path.basename(__file__)),
                    artist='Matplotlib', comment='')
    writer = manimation.writers['ffmpeg'](fps=fps, metadata=metadata)
    writer.setup(fig=fig, dpi=dpi, outfile=savemov)
    for fi in trange(tout.size):  # trange gives us a statusbar if tqdm is installed...
        update(fi)
        writer.grab_frame()
    for _ in range(fps):
        writer.grab_frame()  # 1 second without change (for easier introspection)
    writer.finish()


def average_length_by_weight(Cout, tokg):
    sz = Cout[next(iter(Cout))].size
    nrm = np.zeros(sz)
    val = np.zeros(sz)
    for w, g in tokg['n'].items():
        for k in g:
            val += w*w*Cout[k]
            nrm += w*Cout[k]
    return val/nrm


def total_of(Cout, g):
    if 0 not in g:
        raise ValueError("Probably the wrong token?")
    sz = Cout[next(iter(Cout))].size
    val = np.zeros(sz)
    for num, sks in g.items():
        for sk in sks:
            val += num*Cout[sk]
    return val


def sum_up(Cout, g):
    if 0 in g:
        raise ValueError("Counter token probably the wrong?")
    for k in g:
        if not isinstance(k, int):
            raise ValueError("Wrong group?")
    sz = Cout[next(iter(Cout))].size
    val = np.zeros(sz)
    for num, sks in g.items():
        for sk in sks:
            val += Cout[sk]
    return val
