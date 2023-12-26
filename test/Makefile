.PHONY: all

all: instance bipartite_graph matrix glpk median_heuristic branch_and_cut

clean:
	rm *.exe
	rm *.o

%:
	@if [ $* != "all" ]; then\
		$(CXX) $(CXXFLAGS) -c $*.cpp && \
		$(CXX) $*.o $(LFLAGS) -o $*.exe && \
		./$*.exe; \
	fi