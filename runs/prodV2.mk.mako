<%
import itertools
ter = [-1, -2]
tra = [-2, -3]
concs = [
      {'O2': 0.3e-3, 'N2O': 0.0},
      {'O2': 0.0, 'N2O': 29e-3},
]
all_combos = list(itertools.product(*map(range, map(len, [ter, tra]))))
crn = lambda *args: ('%s-res/crn-%s_{0}_{1}.dat' % (TOKEN, TOKEN)).format(*args)
res = lambda *args: ('%s-res/%s_{0}_{1}_{2}.txt' % (TOKEN, TOKEN)).format(*args)
mp4 = lambda *args: ('%s-res/%s_{0}_{1}_{2}.mp4' % (TOKEN, TOKEN)).format(*args)
tar = lambda *args: ('%s-res/%s_{0}_{1}.tar.xz' % (TOKEN, TOKEN)).format(*args)
%>

MK_${TOKEN}=--n0 ${N0} -n "${N_LEVELS}" -a "[0,1,10,100,1000]" -p "[0,1,10,100,1000]" -o "[0,1,2]" -l "[0,1,10,100,1000,10000]" \
		--reactions rad-rxns.json --substances rad-subst.json
SOLVE_${TOKEN}= --abstol abstol-${TOKEN}.json --settings '{"rtol": 1e-10, "mxsteps": 25000}' --duration 0.0  --varied varied-${TOKEN}.json --npoints 5


.PHONY: all-${TOKEN}

all-${TOKEN}: ${' '.join(['%s %s' % (res(*combo, 0), res(*combo, 1)) for combo in all_combos])}

abstol-${TOKEN}.json: mk_abstol.py
	./mk_abstol.py --outfile $@ --n-levels ${N_LEVELS} --n0 ${N0}

%for (i, j) in all_combos:
${crn(i, j)}: mk_radsys rad-rxns.json rad-subst.json
	mkdir -p ${TOKEN}-res/
	./$< --output $@ $(MK_${TOKEN}) --primary-data '{"loop_scaling_inter": ${ter[i]}, "loop_scaling_intra": ${tra[j]}}'
%for k in range(len(concs)):
${res(i, j, k)}: ${crn(i, j)} solve_ivp varied-${TOKEN}.json abstol-${TOKEN}.json
	./solve_ivp -i $< -o $@ $(SOLVE_${TOKEN}) -c '[{"H2O": 55.4, "H+": 1e-7, "OH-": 1e-7, "n${N0}a0p0o0l0": 1.22e-6, "O2": ${concs[k]['O2']}, "N2O": ${concs[k]['N2O']}}, 0.0]'
${mp4(i, j, k)}: ${res(i, j, k)} anim.py
	./anim.py -r $< --text-conc "N2O,O2,H2O2,HO2,O2-,OH,e-(aq),H" --vmin 1e-15 --vmax 5e-7 --varied varied-${TOKEN}.json --savemov $@

%endfor
${tar(i, j)}: ${' '.join([res(i, j, k) for k in range(len(concs))])}
	tar cJf $@ $^


%endfor
