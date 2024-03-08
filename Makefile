CXX = gcc -std=c++17 -Wall -Werror
CXXFLAGS = -I../src
LFLAGS = -lglpk -lClp -lstdc++ -lm

ifdef PRODUCTION
	CXXFLAGS += -O3 -DNDEBUG
else
	CXXFLAGS += -O0 -g
endif

ifdef VERBOSE
	CXXFLAGS += -DNDEBUG_PRINT
endif

.PHONY: all

all:

clean:
	$(MAKE) -C test clean

test:
	$(MAKE) -C test test.all CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS)" LFLAGS="$(LFLAGS)"

test.%:
	$(MAKE) -C test $* CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS)" LFLAGS="$(LFLAGS)" INPUT=$(INPUT)
