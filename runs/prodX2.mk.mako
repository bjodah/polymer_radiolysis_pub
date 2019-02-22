<%
import itertools
import json
cases = __import__('cases_%s' % TOKEN)
res = lambda **kw: ('%s-res/%s_{mf}_{atmo}_{doserate}.txt' % (TOKEN, TOKEN)).format(**kw)
mp4 = lambda **kw: ('%s-res/%s_{mf}_{atmo}_{doserate}.mp4' % (TOKEN, TOKEN)).format(**kw)

%>

.PHONY: all-${TOKEN}

all-${TOKEN}: ${' '.join([res(mf=mf, atmo=atmo, doserate=doserate) for mf, atmo, doserate in cases.int_combos])}

%for doserate in range(len(cases.doserates)):
varied-${TOKEN}-${doserate}.json: mk_varied.py
	./$< --outfile $@ --tot-dose-Gy 40e3 --f-pulse-Hz 40 --count-wait "[]" --dr-pulse-Gy-s ${cases.doserates[doserate]}
%endfor

%for mf, atmo, doserate in cases.int_combos:
${res(mf=mf, atmo=atmo, doserate=doserate)}: W3-res/crn-W3_0.dat solve_ivp varied-${TOKEN}-${doserate}.json abstol-W3.json
	./solve_ivp -i $< -o $@ --abstol abstol-W3.json --settings '{"rtol": 1e-10, "mxsteps": 25000}' --duration 0.0  --varied varied-${TOKEN}-${doserate}.json  --npoints 5 \
            -c '[{"H2O": 55.4, "H+": 1e-7, "OH-": 1e-7, "n${N0}a0p0o0l0": ${cases.mass_fracs[mf]*998/(410e3)}, "O2": ${cases.atmos[atmo]['O2']}, "N2O": ${cases.atmos[atmo]['N2O']}}, 0.0]'

${mp4(mf=mf, atmo=atmo, doserate=doserate)}: ${res(mf=mf, atmo=atmo, doserate=doserate)} anim.py
	./anim.py -r $< --text-conc "N2O,O2,H2O2,HO2,O2-,OH,e-(aq),H" --vmin 1e-15 --vmax ${10.0 * cases.mass_fracs[mf]*998/(410e3)} --varied varied-${TOKEN}-${doserate}.json --savemov $@
%endfor
