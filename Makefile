CXX ?= g++-8
ifeq ($(CXX), g++)
    COMPILER_ID=GNU
else
    COMPILER_ID=CLANG
endif
LAPACK ?= -llapack -lblas
LDLIBS ?=-lsundials_cvodes -lsundials_nvecserial $(LAPACK) $(EXTRA_LIBS)
# -Wmissing-declarations
WARNINGS ?= \
-Wall \
-Wextra \
-Wredundant-decls \
-Wcast-align \
-Wmissing-include-dirs \
-Wswitch-enum \
-Wswitch-default \
-Winvalid-pch \
-Wredundant-decls \
-Wformat=2 \
-Wmissing-format-attribute \
-Wformat-nonliteral \
-Wodr
INCLUDE ?= -I./src
CXXFLAGS += -Werror -pedantic
CXXFLAGS += -std=c++17 -DKINETGINE_HEADER_ONLY
CXXFLAGS += $(INCLUDE)
CXXFLAGS += $(DEFINES)
CXXFLAGS += $(WARNINGS)
CXXFLAGS += $(EXTRA_FLAGS)
CXXFLAGS += $(shell cmake --find-package -DCMAKE_PREFIX_PATH=$(CMAKE_PREFIX_PATH) -DNAME=SymEngine -DCOMPILER_ID=$(COMPILER_ID) -DLANGUAGE=CXX -DMODE=COMPILE)
LDLIBS += $(shell cmake --find-package -DCMAKE_PREFIX_PATH=$(CMAKE_PREFIX_PATH) -DNAME=SymEngine -DCOMPILER_ID=$(COMPILER_ID) -DLANGUAGE=CXX -DMODE=LINK)

ifeq ($(DEBUG),1)
  CONTEXT=gdb -ex "catch throw" -ex "run" --args
  ifneq (,$(findstring clang++,$(CXX)))  # found clang++
    CXXFLAGS += -g
  else
    CXXFLAGS += -g -ggdb
  endif
else
  CXXFLAGS += -DNDEBUG -O2
endif

TARGETS=mk_odesys solve_ivp rad-res.txt mk_radsys

all: $(TARGETS)

clean:
	rm $(TARGETS)

.PHONY: all clean


sovle_ivp: solve_ivp.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

include rad.mk

mk_radsys: src/mk_radsys.cpp src/ensemble.hpp
	$(CXX) $(CXXFLAGS) -o $@ -DHAVE_KINETGINE $< $(LDLIBS)
