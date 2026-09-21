from __future__ import annotations

import json

import sympy as sp

from bg import (
    BG2CellMLV1,
    BondGraph,
    ComponentType,
    Domain,
    DomainRefiner,
    EquationBuilder,
    SCAPEngine,
    exportBG,
    importBG,
)
from bg.__main__ import build_parser


def make_mass_spring_graph() -> BondGraph:
    """Build the small graph shared by topology, SCAP, and math tests."""
    graph = BondGraph("test_graph")
    graph.add_component("source", type=ComponentType.SE)
    graph.add_component("junction", type=ComponentType.ONE)
    graph.add_component("mass", type=ComponentType.I)
    graph.add_component("spring", type=ComponentType.C)
    graph.add_component("damper", type=ComponentType.R)
    graph.add_bond("source", "junction")
    graph.add_bond("junction", "mass")
    graph.add_bond("junction", "spring")
    graph.add_bond("junction", "damper")
    return graph


def test_graph_editing_and_json_round_trip(tmp_path):
    graph = make_mass_spring_graph()
    assert len(graph.components) == 5
    assert len(graph.bonds) == 4

    graph.delete_bond(graph.components["junction"], graph.components["damper"])
    assert len(graph.bonds) == 3
    graph.add_bond(graph.components["junction"], graph.components["damper"])

    output = tmp_path / "graph.json"
    exportBG(graph, output.name, str(tmp_path))
    restored = importBG(output.name, str(tmp_path))

    assert restored is not None
    assert restored.name == graph.name
    assert set(restored.components) == set(graph.components)
    assert len(restored.bonds) == len(graph.bonds)


def test_domain_refinement_applies_variables_parameters_and_equations():
    graph = BondGraph("refined")
    graph.add_component("spring", type=ComponentType.C)
    refiner = DomainRefiner("domain_catalog.json")

    refiner.refine_graph(
        graph,
        {"spring": (Domain.MECHANICAL_TRANSLATIONAL, "C")},
    )

    spring = graph.components["spring"]
    assert spring.domain is Domain.MECHANICAL_TRANSLATIONAL
    assert spring.ports["1"].effort.physical_quantity is not None
    assert spring.ports["1"].quantity.physical_quantity is not None
    assert spring.constitutive_equations
    assert spring.parameters


def test_scap_assigns_all_power_bond_causalities():
    graph = make_mass_spring_graph()
    system_type = SCAPEngine(graph).run()

    assert system_type.name in {"ODE", "DAE_DERIVATIVE", "DAE_ALGEBRAIC", "DAE_MIXED"}
    for bond in graph.bonds:
        assert bond.validate_causality()


def test_equation_builder_solves_ode_expression():
    graph = BondGraph("equations")
    builder = EquationBuilder(graph)

    equation = builder._process_ode_expression(
        sp.parse_expr("ode(x, t) - a*x")
    )

    assert equation is not None
    assert equation.dependent_symbol == "x"
    assert equation.voi == "t"
    assert equation.infix_rhs == "a*x"


def test_cellml_generation_includes_model_component():
    graph = BondGraph("cellml_model")
    graph.equations = []
    model = BG2CellMLV1(graph)

    assert model.tag == "model"
    component = model.find("component")
    assert component is not None
    assert component.attrib["name"] == "cellml_model"


def test_cli_parser_accepts_refinement_command():
    parser = build_parser()
    args = parser.parse_args(
        [
            "refine",
            "input.json",
            "--component",
            "spring=MECHANICAL_TRANSLATIONAL:C",
            "--output",
            "output.json",
        ]
    )

    assert args.command == "refine"
    assert args.component == [
        ("spring", Domain.MECHANICAL_TRANSLATIONAL, "C")
    ]


def test_exported_json_is_valid(tmp_path):
    output = tmp_path / "empty.json"
    exportBG(BondGraph("empty"), output.name, str(tmp_path))

    with output.open(encoding="utf-8") as file:
        data = json.load(file)

    assert data["name"] == "empty"
    assert data["components"] == []
    assert data["bonds"] == []