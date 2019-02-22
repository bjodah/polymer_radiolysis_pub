rad-crn.dat: mk_odesys rad-rxns.json rad-subst.json
	./$< --reactions rad-rxns.json --substances rad-subst.json --output $@ --verbose --identify-quads >rad-crn.log

rad-res.txt: solve_ivp rad-crn.dat
	./$< -i rad-crn.dat -o $@ -p '{"doserate": 42.0}' \
	-c '[{"H2O": 55.4, "H+": 1e-7, "OH-": 1e-7}, 0.0]' \
	--abstol '[{"H2O": 55e-6, "H+": 1e-13, "OH-": 1e-13, "H2O2": 1e-13, "H2": 1e-13, "O2": 1e-13}, 1e-19]' \
	--settings '{"rtol": 1e-10, "mxsteps": 2000}'

rad-res.png: plot.py rad-res.txt
	MPLBACKEND=Agg ./$< --savefig $@ --results $(word 2,$^)
