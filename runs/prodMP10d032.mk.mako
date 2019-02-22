<%
import itertools
import json
cases = __import__('cases_%s' % TOKEN)
crn = lambda **kw: ('MP10d40-res/crn-MP10d40_{scal}.dat').format(**kw)
res = lambda **kw: ('%s-res/%s_{scal}_{mf}_{atmo}.txt' % (TOKEN, TOKEN)).format(**kw)
mp4 = lambda **kw: ('%s-res/%s_{scal}_{mf}_{atmo}.mp4' % (TOKEN, TOKEN)).format(**kw)

N0=100

%>

SOLVE_ARGS_${TOKEN}= --abstol abstol-MP10d40mf.json --settings '{"rtol": 1e-10, "mxsteps": 25000}' --duration 0.0 --varied varied-${TOKEN}.json --npoints 300


.PHONY: all-${TOKEN}

all-${TOKEN}: ${' '.join([' '.join([res(scal=scal, mf=mf, atmo=atmo) for mf, atmo in cases.int_combos]) for scal, in cases.crn_combos])}

%for scal, in cases.crn_combos:
%for mf, atmo in cases.int_combos:
${res(scal=scal, mf=mf, atmo=atmo)}: ${crn(scal=scal)} solve_ivp varied-${TOKEN}.json abstol-MP10d40mf.json
	mkdir -p ${TOKEN}-res/
	./solve_ivp -i $< -o $@ $(SOLVE_ARGS_${TOKEN}) -c '[{"H2O": 55.4, "H+": 1e-7, "OH-": 1e-7, "n${N0}a0p0o0l0": ${cases.mass_fracs[mf]*998/(float(MW_kDa)*1e3)}, "O2": ${cases.atmos[atmo]['O2']}, "N2O": ${cases.atmos[atmo]['N2O']}}, 0.0]'
${mp4(scal=scal, mf=mf, atmo=atmo)}: ${res(scal=scal, mf=mf, atmo=atmo)} anim.py
	./anim.py -r $< --text-conc "N2O,O2,H2O2,HO2,O2-,OH,e-(aq),H" --vmin 1e-15 --vmax ${10.0 * cases.mass_fracs[mf]*998/(float(MW_kDa)*1e3)} --varied varied-${TOKEN}.json --savemov $@

%endfor
%endfor
