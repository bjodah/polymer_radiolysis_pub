<%
import itertools
mass_fracs = [5e-4, 5e-3]
concs = [
      {'O2': 0.3e-3, 'N2O': 0.0},
      {'O2': 0.0, 'N2O': 29e-3},
]
crn = lambda: ('%s-res/crn-%s.dat' % (TOKEN, TOKEN))
res = lambda *args: ('%s-res/%s_{0}_{1}.txt' % (TOKEN, TOKEN)).format(*args)
mp4 = lambda *args: ('%s-res/%s_{0}_{1}.mp4' % (TOKEN, TOKEN)).format(*args)
%>

MK_${TOKEN}=--n0 ${N0} -n "${N_LEVELS}" -a "[0,1,10,100,1000]" -p "[0,1,10,100,1000]" -o "[0,1,2]" -l "[0,1,10,100,1000,10000]" \
		--reactions rad-rxns.json --substances rad-subst.json
SOLVE_${TOKEN}= --abstol abstol-${TOKEN}.json --settings '{"rtol": 1e-10, "mxsteps": 25000}' --duration 0.0  --varied varied-${TOKEN}.json --npoints 5


.PHONY: all-${TOKEN}

all-${TOKEN}: ${' '.join([res(i, j) for i in range(len(mass_fracs)) for j in range(len(concs))])}

abstol-${TOKEN}.json: mk_abstol.py
	./mk_abstol.py --outfile $@ --n-levels ${N_LEVELS} --n0 ${N0}

${crn()}: mk_radsys rad-rxns.json rad-subst.json
	mkdir -p ${TOKEN}-res/
	./$< --output $@ $(MK_${TOKEN}) --primary-data '{"loop_scaling_inter": -1, "loop_scaling_intra": -2}'
%for i, mf in enumerate(mass_fracs):
%for j in range(len(concs)):
${res(i, j)}: ${crn()} solve_ivp varied-${TOKEN}.json abstol-${TOKEN}.json
	./solve_ivp -i $< -o $@ $(SOLVE_${TOKEN}) -c '[{"H2O": 55.4, "H+": 1e-7, "OH-": 1e-7, "n${N0}a0p0o0l0": ${mf*998/(410e3)}, "O2": ${concs[j]['O2']}, "N2O": ${concs[j]['N2O']}}, 0.0]'
${mp4(i, j)}: ${res(i, j)} anim.py
	./anim.py -r $< --text-conc "N2O,O2,H2O2,HO2,O2-,OH,e-(aq),H" --vmin 1e-15 --vmax ${10.0 * mf*998/(410e3)} --varied varied-${TOKEN}.json --savemov $@

%endfor
%endfor
