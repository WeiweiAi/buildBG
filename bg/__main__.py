"""Command-line entry point for common bond-graph workflows."""

from __future__ import annotations

import argparse
from pathlib import Path

from bg import (
    BondGraph,
    ComponentType,
    ConnectionType,
    Domain,
    DomainRefiner,
    exportBG,
    importBG,
    drawBG,
    print_bond_table,
    SCAPEngine,
    EquationBuilder,
    BG2CellMLV1,
    write_cellmlV1
)



def _component_type(value: str) -> ComponentType | str:
    """Convert a known component type name to its enum, preserving custom types."""
    return ComponentType.__members__.get(value.upper(), value)


def _domain(value: str) -> Domain:
    """Convert a domain name to its enum and reject unknown domains early."""
    try:
        return Domain[value.upper()]
    except KeyError as error:
        valid = ", ".join(domain.name for domain in Domain)
        raise argparse.ArgumentTypeError(
            f"Unknown domain '{value}'. Choose from: {valid}."
        ) from error


def _connection_type(value: str) -> ConnectionType:
    """Convert a connection type name to its enum for graph editing."""
    try:
        return ConnectionType[value.upper()]
    except KeyError as error:
        valid = ", ".join(connection.name for connection in ConnectionType)
        raise argparse.ArgumentTypeError(
            f"Unknown connection type '{value}'. Choose from: {valid}."
        ) from error


def _load_graph(path: Path) -> BondGraph:
    """Load a JSON bond graph and provide a consistent command error."""
    try:
        return importBG(str(path))
    except (OSError, ValueError, KeyError) as error:
        raise RuntimeError(f"Could not load bond graph '{path}': {error}") from error


def _write_graph(graph: BondGraph, path: Path) -> None:
    """Serialize a graph to JSON, creating its parent directory when needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    exportBG(graph, str(path))


def _parse_refinement(value: str) -> tuple[str, Domain, str]:
    """Parse `component=domain[:template]` refinement syntax."""
    try:
        component, specification = value.split("=", 1)
        domain_name, template = (
            specification.split(":", 1)
            if ":" in specification
            else (specification, "")
        )
        if not component or not domain_name:
            raise ValueError
        return component, _domain(domain_name), template
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "Refinement must use COMPONENT=DOMAIN[:TEMPLATE]."
        ) from error


def command_edit(args: argparse.Namespace) -> None:
    """Apply component and bond edits to a graph, then serialize the result."""
    graph = _load_graph(args.input) if args.input else BondGraph(args.name)

    for component_name in args.add_component:
        name, type_name = component_name.split("=", 1)
        graph.add_component(name, type=_component_type(type_name))
    for bond in args.add_bond:
        source, target = bond.split("=", 1)
        graph.add_bond(source, target, type=args.connection_type)
    for bond in args.delete_bond:
        source, target = bond.split("=", 1)
        graph.delete_bond(source, target)

    _write_graph(graph, args.output)


def command_serialize(args: argparse.Namespace) -> None:
    """Re-serialize a graph or convert it to a CellML V1 model."""
    graph = _load_graph(args.input)
    if args.cellml:
        args.cellml.parent.mkdir(parents=True, exist_ok=True)
        write_cellmlV1(BG2CellMLV1(graph), str(args.cellml))
    else:
        _write_graph(graph, args.output)


def command_refine(args: argparse.Namespace) -> None:
    """Apply catalog-backed domain metadata to selected graph components."""
    graph = _load_graph(args.input)
    refinement_map = {
        name: (domain, template or None)
        for name, domain, template in args.component
    }
    DomainRefiner(str(args.catalog)).refine_graph(graph, refinement_map)
    _write_graph(graph, args.output)


def command_scap(args: argparse.Namespace) -> None:
    """Run sequential causality assignment and serialize the classified graph."""
    graph = _load_graph(args.input)
    system_type = SCAPEngine(graph).run()
    _write_graph(graph, args.output)
    print(f"System type: {system_type.name}")


def command_math(args: argparse.Namespace) -> None:
    """Generate, translate, and serialize network or constitutive equations."""
    graph = _load_graph(args.input)
    builder = EquationBuilder(graph)
    equations = []
    if args.kind in ("network", "all"):
        equations.extend(builder.generate_network_equations())
    if args.kind in ("constitutive", "all"):
        equations.extend(builder._generate_constitutive_equations())
    graph.equations = builder.translate_to_symbols(equations)
    _write_graph(graph, args.output)
    print(f"Generated {len(graph.equations)} equation(s).")


def command_visualize(args: argparse.Namespace) -> None:
    """Print a bond table and/or render the graph through Graphviz."""
    graph = _load_graph(args.input)
    if args.table:
        print_bond_table(graph)
    if args.render:
        drawBG(graph, filename=str(args.render), format=args.format, view=args.view)


def build_parser() -> argparse.ArgumentParser:
    """Build the command parser and all high-level workflow subcommands."""
    parser = argparse.ArgumentParser(description="Build and manipulate bond graphs.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    edit = subparsers.add_parser("edit", help="Add or delete components and bonds.")
    edit.add_argument("--input", type=Path, help="Existing JSON graph; omit to start empty.")
    edit.add_argument("--output", type=Path, required=True, help="Output JSON graph.")
    edit.add_argument("--name", default="bond_graph", help="Name for a new graph.")
    edit.add_argument("--add-component", action="append", default=[], metavar="NAME=TYPE")
    edit.add_argument("--add-bond", action="append", default=[], metavar="SOURCE=TARGET")
    edit.add_argument("--delete-bond", action="append", default=[], metavar="SOURCE=TARGET")
    edit.add_argument("--connection-type", type=_connection_type, default=ConnectionType.POWER_BOND)
    edit.set_defaults(func=command_edit)

    serialize = subparsers.add_parser("serialize", help="Convert or re-serialize a graph.")
    serialize.add_argument("input", type=Path, help="Input JSON graph.")
    serialize_outputs = serialize.add_mutually_exclusive_group(required=True)
    serialize_outputs.add_argument("--output", type=Path, help="Output JSON graph.")
    serialize_outputs.add_argument("--cellml", type=Path, help="Write a CellML V1 model instead of JSON.")
    serialize.set_defaults(func=command_serialize)

    refine = subparsers.add_parser("refine", help="Apply domain catalog metadata.")
    refine.add_argument("input", type=Path)
    refine.add_argument("--catalog", type=Path, default=Path("domain_catalog.json"))
    refine.add_argument("--component", action="append", required=True, type=_parse_refinement, metavar="NAME=DOMAIN[:TEMPLATE]")
    refine.add_argument("--output", type=Path, required=True)
    refine.set_defaults(func=command_refine)

    scap = subparsers.add_parser("scap", help="Assign causality and classify the system.")
    scap.add_argument("input", type=Path)
    scap.add_argument("--output", type=Path, required=True)
    scap.set_defaults(func=command_scap)

    math = subparsers.add_parser("math", help="Generate symbolic graph equations.")
    math.add_argument("input", type=Path)
    math.add_argument("--kind", choices=("network", "constitutive", "all"), default="all")
    math.add_argument("--output", type=Path, required=True)
    math.set_defaults(func=command_math)

    visualize = subparsers.add_parser("visualize", help="Print or render a bond graph.")
    visualize.add_argument("input", type=Path)
    visualize.add_argument("--table", action="store_true", help="Print bond causality table.")
    visualize.add_argument("--render", type=Path, help="Graphviz output path without extension.")
    visualize.add_argument("--format", default="png", help="Graphviz output format.")
    visualize.add_argument("--view", action="store_true", help="Open the rendered file.")
    visualize.set_defaults(func=command_visualize)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Parse command-line arguments, execute one workflow, and return an exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        args.func(args)
    except (RuntimeError, ValueError, OSError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
