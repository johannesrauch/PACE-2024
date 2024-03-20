# PACE
## Dependencies
- Modern version of `gcc`; we use `c++17`.
- `make`.
- [glpk](https://www.gnu.org/software/glpk/). `make` expects the glpk library in the standard folder `/usr/local/include`.

## Compiling
### Tests
The command
```
make test.bipartite_graph
```
runs the bipartite_graph.hpp test, for instance.
Run 
```
make test.all
```
for all tests (except the solver tests with large instances).
The command
```
make clean
```
deletes all executables in the test folder.

## Formatting
Please use clang-format with
```
{ BasedOnStyle: Google, IndentWidth: 4, ColumnLimit: 0 }
```
options.

Please use 
```
#ifndef PACE_[YOUR FILE NAME IN CAPS WITH UNDERSCORES]_HPP
#define PACE_[YOUR FILE NAME IN CAPS WITH UNDERSCORES]_HPP

// your code here

#endif
```
to avoid multiple header inclusions.
Please use the namespace `pace` for the project.

## Add files
```
cd existing_repo
git remote add origin https://gitlab.uni-ulm.de/zhh72/pace.git
git branch -M main
git push -uf origin main
```

## Roadmap
- [ ] Switch linear program solver to [CLP](https://github.com/coin-or/Clp).
- [ ] Implement improvements (todos) mentioned in `branch_and_cut::solve`.

## Authors
Johannes Rauch.

## License
[MIT License](https://en.wikipedia.org/wiki/MIT_License).

## Project status
Work in progress.
