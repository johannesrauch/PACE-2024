.PHONY: all

all: input bipartite_graph matrix crossing_number glpk median_heuristic topological_sort

clean:
	rm *.exe
	rm *.o

%:
	@if [ $* != "all" ]; then\
		$(CXX) $(CXXFLAGS) -c $*.cpp && \
		$(CXX) $*.o $(LFLAGS) -o $*.exe && \
		./$*.exe; \
	fi