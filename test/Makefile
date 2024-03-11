.PHONY: all

all: input bipartite_graph matrix crossing_number glpk median_heuristic topological_sort branch_and_cut

clean:
	rm *.exe
	rm *.o
	rm *.log

%:
	@if [ $* != "all" ]; then \
		$(CXX) $(CXXFLAGS) -c $*.cpp && \
		$(CXX) $*.o $(LFLAGS) -o $*.exe && \
		./$*.exe $(INPUT) > $*.log; \
	fi