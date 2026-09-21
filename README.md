# buildBG

`buildBG` is a Python toolkit for building, refining, analyzing, serializing, and visualizing bond-graph models.

## Installation

Requires Python 3.12 or later. With [uv](https://docs.astral.sh/uv/):

```powershell
uv sync
```

The command-line interface is available as:

```powershell
uv run bg --help
```

## Commands

```text
edit       Add or delete components and bonds
serialize  Convert or re-serialize a graph, including CellML output
refine     Apply physical-domain metadata from domain_catalog.json
scap       Assign causality and classify the system
math       Generate symbolic network or constitutive equations
visualize  Print a bond table or render a Graphviz image
```

Examples:

```powershell
# Create a graph from components and a bond
uv run bg edit --name example --output example.json `
	--add-component Source=SE `
	--add-component Junction=ONE `
	--add-bond Source=Junction

# Apply domain metadata
uv run bg refine example.json `
	--component Source=MECHANICAL_TRANSLATIONAL:Se `
	--output refined.json

# Assign causality and generate equations
uv run bg scap refined.json --output causal.json
uv run bg math causal.json --kind all --output equations.json

# Convert a graph to CellML
uv run bg serialize equations.json --cellml model.cellml

# Print causality or render with Graphviz
uv run bg visualize causal.json --table
uv run bg visualize causal.json --render bond_graph
```

## Python API

The package is importable as `bg`:

```python
from bg import BondGraph, ComponentType

graph = BondGraph("example")
graph.add_component("source", type=ComponentType.SE)
```

Domain templates are defined in `bg/domain_catalog.json`; example models and data are in `data/`.
