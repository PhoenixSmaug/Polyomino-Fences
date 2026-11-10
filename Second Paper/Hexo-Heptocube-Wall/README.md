The solvers use NumPy for grid data structures and polycubetools, a framework for parallel greedy search.

Create a virtual environment and install the required packages into it:

```sh
python3 -m venv venv
venv/bin/python -m pip install numpy polycubetools
```

The virtual environment keeps these dependencies separate from your system-wide Python installation. There is no need to activate it when invoking Python through `venv/bin/python`.

Run the hexocube solver with:

```sh
venv/bin/python hexo_greedy.py
```
