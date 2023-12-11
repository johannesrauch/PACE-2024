all: instance bipartite_graph crossing_number_matrix

%:
	@if [ $* != "all" ]; then\
		$(CXX) $(CXXFLAGS) -o $*.exe $*.cpp; \
		./$*.exe; \
	fi