<%
import itertools
import json
cases = __import__('cases_%s' % TOKEN)
res = lambda **kw: ('%s-res/%s_{mf}_{atmo}_{freq}.txt' % (TOKEN, TOKEN)).format(**kw)
mp4 = lambda **kw: ('%s-res/%s_{mf}_{atmo}_{freq}.mp4' % (TOKEN, TOKEN)).format(**kw)

%>

.PHONY: all-${TOKEN}

all-${TOKEN}: ${' '.join([res(mf=mf, atmo=atmo, freq=freq) for mf, atmo, freq in cases.int_combos])}

%for freq in range(len(cases.freqs)):
varied-${TOKEN}-${freq}.json: mk_varied.py
	./$< --outfile $@ --tot-dose-Gy 40e3 --f-pulse-Hz ${cases.freqs[freq]} --count-wait "[]"
%endfor

%for mf, atmo, freq in cases.int_combos:
${res(mf=mf, atmo=atmo, freq=freq)}: W3-res/crn-W3_0.dat solve_ivp varied-${TOKEN}-${freq}.json abstol-W3.json
	./solve_ivp -i $< -o $@ --abstol abstol-W3.json --settings '{"rtol": 1e-10, "mxsteps": 25000}' --duration 0.0  --varied varied-${TOKEN}-${freq}.json  --npoints 5 \
            -c '[{"H2O": 55.4, "H+": 1e-7, "OH-": 1e-7, "n${N0}a0p0o0l0": ${cases.mass_fracs[mf]*998/(410e3)}, "O2": ${cases.atmos[atmo]['O2']}, "N2O": ${cases.atmos[atmo]['N2O']}}, 0.0]'

${mp4(mf=mf, atmo=atmo, freq=freq)}: ${res(mf=mf, atmo=atmo, freq=freq)} anim.py
	./anim.py -r $< --text-conc "N2O,O2,H2O2,HO2,O2-,OH,e-(aq),H" --vmin 1e-15 --vmax ${10.0 * cases.mass_fracs[mf]*998/(410e3)} --varied varied-${TOKEN}-${freq}.json --savemov $@
%endfor
