.PHONY: all
all: ensemble.hpp # unsemble.cpp

%: render.py %.mako
	python3 $< --target $@
