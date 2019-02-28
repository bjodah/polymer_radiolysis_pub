.PHONY: all
all: mk_radsys solve_ivp
CXXFLAGS=-O2 -std=c++17 -DKINETGINE_HEADER_ONLY -L. $(shell $${PYTHON:-python3} -c "from pycvodes import get_include as gi; print(' '.join('-I%s' % I for I in [gi(), '.']))")
LDLIBS=-lsymengine $(shell $${PYTHON:-python3} -c "from pycvodes import config as c; print(' '.join('-l%s' % l for l in c['SUNDIALS_LIBS'].split(',') + c['LAPACK'].split(',')))")

%: render.py %.mako
	python3 $< --target $@

mk_radsys: mk_radsys.cpp ensemble.hpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDLIBS)

solve_ivp: solve_ivp.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDLIBS)
