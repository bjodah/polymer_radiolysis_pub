<%
import itertools
import json
cases = __import__('cases_%s' % TOKEN)
crn = lambda **kw: ('%s-res/crn-%s_{scal}.dat' % (TOKEN, TOKEN)).format(**kw)
res = lambda **kw: ('%s-res/%s_{scal}_{mf}_{atmo}.txt' % (TOKEN, TOKEN)).format(**kw)
mp4 = lambda **kw: ('%s-res/%s_{scal}_{mf}_{atmo}.mp4' % (TOKEN, TOKEN)).format(**kw)

%>

MK_${TOKEN}=--n0 ${N0} -n "${N_LEVELS}" -a "[0,1,10,100,1000]" -p "[0,1,10,100,1000]" -o "[0,1,2]" -l "[0,1,10,100,1000,10000]" \
		--reactions rad-rxns.json --substances rad-subst.json
SOLVE_${TOKEN}= --abstol abstol-${TOKEN}.json --settings '{"rtol": 1e-10, "mxsteps": 25000}' --duration 0.0  --varied varied-${TOKEN}.json --npoints 5



.PHONY: all-${TOKEN}

all-${TOKEN}: ${' '.join([' '.join([res(scal=scal, mf=mf, atmo=atmo) for mf, atmo in cases.int_combos]) for scal, in cases.crn_combos])}

abstol-${TOKEN}.json: mk_abstol.py
	./mk_abstol.py --outfile $@ --n-levels ${N_LEVELS} --n0 ${N0}

%for scal, in cases.crn_combos:
${crn(scal=scal)}: mk_radsys rad-rxns.json rad-subst.json
	mkdir -p ${TOKEN}-res/
	./$< --output $@ $(MK_${TOKEN}) --primary-data '${json.dumps(cases.scals[scal])}'
%for mf, atmo in cases.int_combos:
${res(scal=scal, mf=mf, atmo=atmo)}: ${crn(scal=scal)} solve_ivp varied-${TOKEN}.json abstol-${TOKEN}.json
	./solve_ivp -i $< -o $@ $(SOLVE_${TOKEN}) -c '[{"H2O": 55.4, "H+": 1e-7, "OH-": 1e-7, "n${N0}a0p0o0l0": ${cases.mass_fracs[mf]*998/(410e3)}, "O2": ${cases.atmos[atmo]['O2']}, "N2O": ${cases.atmos[atmo]['N2O']}}, 0.0]'
${mp4(scal=scal, mf=mf, atmo=atmo)}: ${res(scal=scal, mf=mf, atmo=atmo)} anim.py
	./anim.py -r $< --text-conc "N2O,O2,H2O2,HO2,O2-,OH,e-(aq),H" --vmin 1e-15 --vmax ${10.0 * cases.mass_fracs[mf]*998/(410e3)} --varied varied-${TOKEN}.json --savemov $@

%endfor
%endfor
