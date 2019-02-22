parameter = dict(XB1='freqs', XB2='doserates_Gy_s', XB3='tpulses_s')
offset = dict(XB1=0, XB2=1, XB3=0)[TOKEN]
param_pretty = dict(freqs='Frequency', doserates_Gy_s='Doserate', tpulses_s='Pulse length')
param_unit = dict(freqs='Hz', doserates_Gy_s='Gy/s', tpulses_s='s')
vark = parameter[TOKEN]
varl = dict(XB1=[], XB2=[32e6], XB3=[5e-6])[TOKEN] + getattr(cases, vark)
from chempy.printing import number_to_scientific_latex

for mfi in range(len(cases.mass_fracs)):
    for xi in range(len(varl)-offset):
        cx, cy, ip1d = interpols[xi]
        ppm = cases.mass_fracs[mfi]*1e6
        # assert int(ppm) == 5000
        air, nit = all_air[mfi, xi], all_nit[mfi, xi]
        t_end_s = air['tout'][-1]
        print('t_end = %.5f s' % t_end_s)
        fig, ax = plt.subplots(1, 1)
        final_length_by_weight={}
        for lbl, cont in list(zip(['air', 'nitrous oxide'], [air, nit])):
            n0 = int(cont['sk0'][1:].split('a')[0])  # intial lenght
            ax.axhline(n0, ls='--', lw=0.5, alpha=0.5, c='k')
            avglen = average_length_by_weight(cont['Cout'], cont['tokg'])
            ax.plot(ip1d(cont['tout'])/1e3, avglen, label=lbl)
            final_length_by_weight[lbl] = avglen[-1]
        ax.set_ylabel('Chain length / number of segments')
        ax.set_xlabel('Dose / kGy')
        ax.legend()
        ax.set_ylim([90, None])

        varv = varl[xi]
        ax.set_title('{} = ${}$ {}'.format(param_pretty[vark], number_to_scientific_latex(varv), param_unit[vark]))
        fnam = lambda k: 'poly_%s_%s/{}ppm_{}_{}'.format(ppm, vark, varv) % (TOKEN, k.replace('(', '').replace(')', ''))
        def savefig(figur, key):
            path = '%s.png' % fnam(key)
            if not os.path.exists(os.path.dirname(path)):
                os.mkdir(os.path.dirname(path))
            figur.savefig(path, dpi=300)

        savefig(fig, 'size')
        dur, varied = json.load(open(_varied(xi)))
        n_pulse = len(dur)/2
        dr = varied['doserate'][0]
        with open('%s_%dppm.txt' % (fnam('size'), ppm), 'wt') as fh:
            fh.write("""t_end_s={}, freq_Hz={}, t_pulse_s={}, t_wait_s={}, doserate_Gy_s={:.4g},
final_length_by_weight_air={}, final_length_by_weight_N2O={}""".format(
                t_end_s, n_pulse/t_end_s, dur[0], dur[1], dr, final_length_by_weight['air'], final_length_by_weight['nitrous oxide']
            ))

        for specie, ylim in [('H2O2', [0, 2e-3]), ('e-(aq)', [0, 2e-5])]:
            fig2, ax2 = plt.subplots(1, 1)
            for lbl, cont in list(zip(['air', 'nitrous oxide'], [air, nit])):
                ax2.plot(ip1d(cont['tout'])/1e3, cont['Cout'][specie], label=lbl, alpha=0.7)
            ax2.set_ylabel(r'[$\mathrm{%s}$] / M' % {'H2O2': 'H_2O_2', 'e-(aq)': 'e^-(aq)'}[specie])
            ax2.set_xlabel('Dose / kGy')
            ax2.set_ylim(ylim)
            ax2.legend()
            savefig(fig2, specie)

        # unsaturated
        fig3, ax3 = plt.subplots(1, 1)
        for lbl, cont in list(zip(['air', 'nitrous oxide'], [air, nit])):
            ax3.plot(ip1d(cont['tout']), sum_up(cont['Cout'], cont['ctrg']['unsaturated'])/cont['ic0'], label=lbl)
        ax3.set_ylabel('Number of double bonds per polymer')
        ax3.set_xlabel('Dose / kGy')
        ax3.set_ylim([0, 400])
        ax3.legend()
        savefig(fig3, 'unsaturated')
