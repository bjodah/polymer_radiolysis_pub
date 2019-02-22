#!/usr/bin/env python3


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

def main(radsys_factory='mk_radsys', MW_kDa=410.0, mass_fractions="[1e-9,5e-4,1e-3,2e-3,5e-3]", token='T', append=False,
         npoints=2, n0=10, n='[1,10,100]', a='[0,1,10,100]', p='[0,1,10,100]', o='[0,1,2]', l='[0,1,10,100,1000]', primary_data=""):
    if not radsys_factory:
        raise ValueError("You must specify radsys_factory, e.g. 'mk_radsys_mod'")
    concs = []
    for mf in eval(mass_fractions):
        poly_conc = mf*998/(MW_kDa*1e3)
        concs.append(dict(N2O=25e-3, O2=0, polymer=poly_conc))
        concs.append(dict(N2O=0e-3, O2=0.26e-3, polymer=poly_conc))

    fname = 'prod{}.mk'.format(token)
    variables = dict(
        PROD_TOK=token,
        RADSYS_FACTORY=radsys_factory,
        concenrations=concs,
        N0=n0,
        LEVELS_N=n,
        LEVELS_A=a,
        LEVELS_P=p,
        LEVELS_O=o,
        LEVELS_L=l,
        NPOINTS=npoints,
        PRIMARY_DATA=primary_data
    )
    open(fname, 'wt').write(render_mako_template(open('prodT.mk.mako'), variables).decode('utf-8'))
    if append:
        open('Makefile', 'a').write('include {}\n'.format(fname))


if __name__ == '__main__':
    argh.dispatch_command(main)
