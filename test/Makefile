.PHONY: all

all: input bipartite_graph matrix crossing_number glpk clp median_heuristic topological_sort branch_and_cut

clean:
	rm *.exe
	rm *.o

%:
	@if [ $* != "all" ]; then\
		$(CXX) $(CXXFLAGS) -c $*.cpp && \
		$(CXX) $*.o $(LFLAGS) -o $*.exe && \
		./$*.exe; \
	fi