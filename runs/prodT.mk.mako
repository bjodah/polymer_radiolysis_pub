<%doc>
This is a template. Render using ``render.py``.
</%doc>

prod${PROD_TOK}-crn.dat: ${RADSYS_FACTORY} rad-rxns.json rad-subst.json
	./$< --n0 ${N0} -n '${LEVELS_N}' -a '${LEVELS_A}' -p '${LEVELS_P}' -o '${LEVELS_O}' -l '${LEVELS_L}' \
	--reactions rad-rxns.json --substances rad-subst.json --output $@ \
        --primary-data '${PRIMARY_DATA}'

%for idx, conc in enumerate(concenrations):
prod${PROD_TOK}-${idx}.txt: prod${PROD_TOK}-crn.dat varied-${PROD_TOK}.json
	./solve_ivp -i $< -o $@ --verbose --npoints ${NPOINTS} \
	-c '[{"H2O": 55.4, "H+": 1e-7, "OH-": 1e-7, "n${N0}a0p0o0l0": ${conc["polymer"]}, "O2": ${conc["O2"]}, "N2O": ${conc["N2O"]}}, 0.0]' \
	--abstol '[{"H2O": 55e-6, "H+": 1e-13, "OH-": 1e-13, "H2O2": 1e-13, "H2": 1e-13, "O2": 1e-13}, 1e-19]' \
	--settings '{"rtol": 1e-10, "mxsteps": 25000}' --duration 0.0  --varied varied-${PROD_TOK}.json

prod${PROD_TOK}-${idx}.png: plot.py prod${PROD_TOK}-${idx}.txt
	./$< --savefig $@ --results $(word 2,$^)
%endfor
