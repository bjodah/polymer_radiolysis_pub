rad-crn.dat: rad-rxns.json rad-subst.json
	mk_radsys --reactions rad-rxns.json --substances rad-subst.json --n0 2 --output $@ --verbose --identify-quads >rad-crn.log

rad-res.txt: rad-crn.dat
	solve_ivp -i rad-crn.dat -o $@ -p '{"doserate": 42.0}' \
	-c '[{"H2O": 55.4, "H+": 1e-7, "OH-": 1e-7, "n2a0p0o0l0": 0.01}, 0.0]' \
	--abstol '[{"H2O": 55e-6, "H+": 1e-13, "OH-": 1e-13, "H2O2": 1e-13, "H2": 1e-13, "O2": 1e-13}, 1e-19]' \
	--settings '{"rtol": 1e-10, "mxsteps": 2000}'

rad-res.png: plot.py rad-res.txt
	MPLBACKEND=Agg python3 ./$< --savefig $@ --results $(word 2,$^)
