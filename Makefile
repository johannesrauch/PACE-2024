CXX = g++ -std=c++17 -Wall
CXXFLAGS = -I ../src

ifdef PRODUCTION
	CXXFLAGS += -O3
else
	CXXFLAGS += -g
endif

.PHONY: all

all:

clean:
	rm test/*.exe

test.%:
	$(MAKE) -C test $* CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS)"
