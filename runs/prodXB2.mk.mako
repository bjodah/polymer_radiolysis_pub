<%
import itertools
import json
cases = __import__('cases_%s' % TOKEN)
crn = lambda **kw: ('%s-res/crn-%s_{icrn}.dat' % (TOKEN, TOKEN)).format(**kw)
res = lambda **kw: ('%s-res/%s_{icrn}_{iint}.txt' % (TOKEN, TOKEN)).format(**kw)
mp4 = lambda **kw: ('%s-res/%s_{icrn}_{iint}.mp4' % (TOKEN, TOKEN)).format(**kw)

N0=100
N_LEVELS="[10,100,200,600,2400]"

alkyl_inter_M_s=float(alkyl_inter_M_s_410kDa) * (float(MW_kDa)/410.0)**-0.6
peroxy_inter_M_s=alkyl_inter_M_s

%>

MK_${TOKEN}=--n0 ${N0} -n "${N_LEVELS}" -a "[0,1,2,6,24,120]" -p "[0,1,2,4,8]" -o "[0,1,2]" -l "[0,1,10,100,1000,10000]" \
		--reactions rad-rxns.json --substances rad-subst.json
SOLVE_ARGS_${TOKEN}= --abstol abstol-${TOKEN}.json --settings '{"rtol": 1e-10, "mxsteps": 25000}' --duration 0.0 --npoints 3

abstol-${TOKEN}.json: mk_abstol.py
	./mk_abstol.py --outfile $@ --n-levels ${N_LEVELS} --n0 ${N0}

.PHONY: all-${TOKEN}

all-${TOKEN}: ${' '.join([' '.join([res(icrn=icrn, iint=iint) for iint in range(len(cases.int_indices))]) for icrn in range(len(cases.crn_indices))])}

%for idrat in range(len(cases.doserates_Gy_s)):
varied-${TOKEN}-${idrat}.json: mk_varied.py
	./mk_varied.py --outfile $@   --tot-dose-Gy 40e3 --count-wait "[]" --dr-pulse-Gy-s ${cases.doserates_Gy_s[idrat]} --t-pulse-s 5e-6 --f-pulse-Hz 400.0
%endfor


%for icrn, (iprim,) in enumerate(cases.crn_indices):
${crn(icrn=icrn)}: mk_radsys rad-rxns.json rad-subst.json
	mkdir -p ${TOKEN}-res/
	./$< --output $@ $(MK_${TOKEN}) --primary-data '${json.dumps(dict(cases.prim_dat[iprim], MW_kDa=float(MW_kDa), alkyl_inter_M_s=alkyl_inter_M_s, peroxy_inter_M_s=peroxy_inter_M_s, t_half_loop_closure_intra_s_94kDa=float(t_half_loop_closure_intra_s_94kDa)))}'
%for iint, (imf, iatmo, idrat) in enumerate(cases.int_indices):
${res(icrn=icrn, iint=iint)}: ${crn(icrn=icrn)} solve_ivp varied-${TOKEN}-${idrat}.json abstol-${TOKEN}.json
	./solve_ivp -i $< -o $@ $(SOLVE_ARGS_${TOKEN}) --varied varied-${TOKEN}-${idrat}.json \
-c '[{"H2O": 55.4, "H+": 1e-7, "OH-": 1e-7, "n${N0}a0p0o0l0": ${cases.mass_fracs[imf]*998/(float(MW_kDa)*1e3)}, "O2": ${cases.atmos[iatmo]['O2']}, "N2O": ${cases.atmos[iatmo]['N2O']}}, 0.0]'
${mp4(icrn=icrn, iint=iint)}: ${res(icrn=icrn, iint=iint)} anim.py
	./anim.py -r $< --text-conc ${("N2O," if cases.atmos[iatmo]['N2O'] else "") + "O2,H2O2,HO2,O2-,OH,e-(aq),H"} \
--vmin 1e-15 --vmax ${10.0 * cases.mass_fracs[imf]*998/(float(MW_kDa)*1e3)} --varied varied-${TOKEN}-${idrat}.json --savemov $@

%endfor
%endfor
