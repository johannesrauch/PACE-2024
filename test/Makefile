.PHONY: all

all: bipartite_graph matrix crossings heuristics topological_sort branch_and_cut_tiny

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