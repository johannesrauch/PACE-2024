.PHONY: all

all: instance bipartite_graph crossing_number_matrix glpk

clean:
	rm *.exe
	rm *.o

%:
	@if [ $* != "all" ]; then\
		$(CXX) $(CXXFLAGS) -c $*.cpp; \
		$(CXX) $*.o $(LFLAGS) -o $*.exe; \
		./$*.exe; \
	fi