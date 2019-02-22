<%
import itertools
import json
cases = __import__('cases_%s' % TOKEN)
crn = lambda: 'crn-%s.dat' % TOKEN
res = lambda **kw: ('%s-res/%s_{mf}_{tp}.txt' % (TOKEN, TOKEN)).format(**kw)

# N0=100
# N_LEVELS="[1,10,100,1000,10000]"

N0=100
N_LEVELS="[10,100,1000]"

%>

MK_${TOKEN}=--n0 ${N0} -n "${N_LEVELS}" -a "[0,1,2,3,4,8,16,32,64]" -p "[0,1,10,100]" -o "[0,1,2]" -l "[0,1,10,100]" \
		--reactions rad-rxns.json --substances rad-subst.json --primary-data '{"MW_kDa": ${float(MW_kDa)} ${prim}}'
SOLVE_ARGS_${TOKEN}= --abstol abstol-${TOKEN}.json --settings '{"rtol": 1e-10, "mxsteps": 25000}'


.PHONY: all-${TOKEN}

all-${TOKEN}: ${' '.join([' '.join([res(mf=mf, tp=tp) for mf in range(len(cases.mass_fracs))]) for tp in range(len(cases.pulse_durs))])}

abstol-${TOKEN}.json: mk_abstol.py
	./mk_abstol.py --outfile $@ --n-levels ${N_LEVELS} --n0 ${N0}

${crn()}: mk_radsys rad-rxns.json rad-subst.json
	mkdir -p ${TOKEN}-res/
	./$< --output $@ $(MK_${TOKEN})

%for tp in range(len(cases.pulse_durs)):



%for mf in range(len(cases.mass_fracs)):
${res(mf=mf, tp=tp)}: ${crn()} solve_ivp abstol-${TOKEN}.json
	./solve_ivp -i $< -o $@.pulse $(SOLVE_ARGS_${TOKEN}) --duration ${cases.pulse_durs[tp]} -c '[{"H2O": 55.4, "H+": 1e-7, "OH-": 1e-7, "n${N0}a0p0o0l0": ${cases.mass_fracs[mf]*998/(float(MW_kDa)*1e3)}, "O2": 0.0, "N2O": 29e-3}, 0.0]' --parameters '{"doserate": ${DOSERATE}}'
	./final_conc.py $@.pulse >$@.pulse.final
	./solve_ivp -i $< -o $@.post $(SOLVE_ARGS_${TOKEN}) --duration 300.0 --initial-conditions-path $@.pulse.final --parameters '{"doserate": 0.0}'
	./concat_results.py $@.pulse $@.post >$@

%endfor
%endfor
