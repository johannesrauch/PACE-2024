CXX = gcc -std=c++17 -Wall -Werror
CXXFLAGS = -I../src -I/usr/local/include/highs
LFLAGS = -lhighs -lstdc++ -lm

ifdef PRODUCTION
	CXXFLAGS += -O3 -DNDEBUG
else
	CXXFLAGS += -O0 -g
endif

ifdef VERBOSE
	CXXFLAGS += -DNDEBUG_PRINT
endif

ifdef FAST
	CXXFLAGS += -DPACE_CONST_PROBMEDIAN_ITERATIONS=4
endif

.PHONY: all

all:

clean:
	$(MAKE) -C test clean

test:
	$(MAKE) -C test test.all CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS)" LFLAGS="$(LFLAGS)"

test.%:
	$(MAKE) -C test $* CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS)" LFLAGS="$(LFLAGS)" INPUT=$(INPUT)
