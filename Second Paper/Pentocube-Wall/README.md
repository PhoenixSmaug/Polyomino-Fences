The `dlx_solver.py` uses SSXCC, Donald Knuth's dancing-links solver. Build the executable before running the search. The Makefile target `ssxcc` compiles `src/ssxcc.c` together with `src/gb_flip.c`; `src/gb_flip.c` and `src/gb_flip.h` provide the random-number generator routines required by SSXCC.

On the top level the solver generates pentacube placements and encodes the problem as an exact-cover instance, which is then solved by SSXCC. When a solution is found, `postprocessor.py` translates the selected placements back into full pentacube coordinates and writes the grid JSON.

From the project directory, create the Python virtual environment and install the required packages, activate it, build SSXCC, and run the wrapper:

```sh
make venv
source venv/bin/activate
make ssxcc
python dlx_solver.py
```

The run creates `output/solfound.dlx` with the solver instance and solution, and `output/solfound.json` with the piece coordinates. Building requires `make` and a C compiler; the generated `ssxcc` executable is compiled for the machine where the build runs.
