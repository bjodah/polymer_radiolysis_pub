#!/usr/bin/env python3


import os
import io
import argh
from mako.template import Template
from mako.exceptions import text_error_template


def render_mako_template(source_template, variables):
    try:
        template_str = source_template.read()
    except:
        template_str = io.open(source_template, 'rt', encoding='utf-8').read()

    try:
        rendered = Template(template_str, input_encoding='utf-8',
                            output_encoding='utf-8').render(**variables)
    except:
        print(text_error_template().render())
        raise
    return rendered

def main(
        target='ensemble.hpp',
        k_oxyl_H_intra_s=1.5e5,             # "Jonsson" (1000x 150/s, där 150/s är från jeszka_monte-carlo_2006, p. 867)
        k_oxyl_frag_s=1e7,                  # Scheme 3, p 1804 in Salamone & Bietti, Reaction Pathways, Synlett 2014, 15,
        k_alkyl_O2_M_s=1e9,                 # typical;                              1803-1816, doi: 10.1055/s-033-1341280
        k_OH_1kDa_M_s=1.3e9,                # Fig 3, p 478, doi:10.1002/kin.20575, Bartoszek et al.
        ratio_OH_H=102,                     # k(OH)/k(H) on c-hexane = 102
        MW_kDa=410.0,                       # Sweden 410, Poland 640, Ulanski 1300, Kadublowski 2000, cf. 91 kDa
        MW_monomer_g_mol=111.14,            # N-vinylpyrrolidone
        t_half_loop_closure_intra_s_94kDa=0.008,    # extrapolation Figure 17, p 867, jeszka_monte-carlo_2006
        alkyl_inter_M_s=2.5e7,              # Fig 5, p 182, borgwardt_pulsradiolytische_1969
        peroxy_inter_M_s=2.5e7,             # based on alykl_inter
        N_scal_intra=-2.0625,               # jeszka_monte-carlo_2006, fig 8 (-np.mean([2.3,2.15,2.1,1.7]))
        w_scal_OH=-0.3,                     # Fig 3, p 478, Bartoszek et al. doi:10.1002/kin.20575
        radrad_scal_intra=2.719,          # jeszka_monte-carlo_2006, p. 863, fig 8:
                                             # t½(a=2) = 2e5, t½(a=16)=700 =>
                                             # dlnt½/da = (log(2e5) - log(700))/(log(2) - log(16)) = -2.71948
        disprop_frac=0.38,                  # disprop=0.38 from 0.63/(1 + 0.63) where 0.63 == kd/kc from
                                            # p. 150 General Aspects of the Chemistry of Radicals
        loop_scaling_intra=-0.5,
        loop_scaling_inter=-0.25            # p. 183, borgwardt_pulsradiolytische_1969, U. BORGWARDT, W. SCHNABEL und A. HENGLEIN, Die Makrornolekulare Chemie 127 (1969) 176-184 (Nr. 3134)
):
    _exclude = ['target']
    primary_data = dict([(k, v) for k, v in locals().items() if k not in _exclude and k != '_exclude'])
    variables = dict(
        PRIMARY_DATA=primary_data
    )
    template = target + '.mako'
    if not os.path.exists(template):
        raise ValueError("Template does not exist: %s" % template)
    txt = render_mako_template(open(template), variables).decode('utf-8')
    open(target, 'wt').write(txt)


if __name__ == '__main__':
    argh.dispatch_command(main)
