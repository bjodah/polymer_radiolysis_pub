.PHONY: all
all: mk_radsys solve_ivp
CXXFLAGS=-O2 -std=c++17 -DKINETGINE_HEADER_ONLY
LDFLAGS=-lsymengine -lsundials_cvodes -lsundials_nvecserial -lsundials_sunlinsollapackdense -lsundials_sunlinsollapackband

%: render.py %.mako
	python3 $< --target $@

mk_radsys: mk_radsys.cpp ensemble.hpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

solve_ivp: solve_ivp.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)
