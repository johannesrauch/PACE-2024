# deprecated!

CXX = gcc -std=c++17 -Wall -Werror
CXXFLAGS = -I../src -I../highs/include/highs
LFLAGS = -lhighs -lstdc++ -lm

ifdef RELEASE
	CXXFLAGS += -O3 -DNDEBUG
else
	CXXFLAGS += -O0 -g
endif

ifdef VERBOSE
	CXXFLAGS += -DPACE_DEBUG_PRINT
endif

.PHONY: all

all:

clean:
	$(MAKE) -C test clean

test:
	$(MAKE) -C test test.all CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS)" LFLAGS="$(LFLAGS)"

test.%:
	$(MAKE) -C test $* CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS)" LFLAGS="$(LFLAGS)" INPUT=$(INPUT)
