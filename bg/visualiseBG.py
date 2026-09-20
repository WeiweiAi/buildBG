import graphviz
from .defineBG import BondGraph, ComponentType, JUNCTIONS, ConnectionType

# Console and Graphviz renderers for inspecting a bond graph after construction.
def print_bond_table(bg: BondGraph) -> None:
    """Prints each bond's endpoints and causal direction for quick inspection."""
    print(f"\n--- Causality Summary: {bg.name} ---")
    print(f"{'Bond':<8} | {'Source':<12} | {'Target':<15} | {'Causality'}")
    print("-" * 60)

    for i, b in enumerate(bg.bonds, 1):
        src = b.source.component.name
        tgt = b.target.component.name

        if b.target.causality == True:
            direction = f"{src} |-----> {tgt}"
        elif b.source.causality == True:
            direction = f"{src} <-----| {tgt}"
        else:
            direction = f"{src} ------- {tgt} (UNASSIGNED)"

        print(f"Bond {i:<3} | {src:<12} | {tgt:<15} | {direction}")

def drawBG(bg: BondGraph, filename: str = "bond_graph", format: str = "png", view: bool = True) -> graphviz.Digraph:
        """
    Renders the bond graph from source to target with formal causal strokes.
        - Power flow arrow points from source to target.
        - Causal stroke is drawn at the effort-receiving end.

    The rendered file is written through Graphviz; the returned graph is
    also available to callers that want to inspect or render it again.
        """
        dot = graphviz.Digraph(name=filename, comment="Bond Graph Visualization")
        dot.attr(rankdir="LR", nodesep="0.6", ranksep="0.8")

        # Clean textbook node styling (no bounding boxes)
        dot.attr("node", shape="plaintext", fontname="Helvetica-Bold", fontsize="14")

        # Render Component Nodes
        for comp_name, comp in bg.components.items():
            label=''
            if comp.type in JUNCTIONS:
                if comp.type == ComponentType.ONE:
                    label = "1"
                elif comp.type == ComponentType.ZERO:
                    label = "0"
                elif comp.type == ComponentType.XONE:
                    label = "X1"
                elif comp.type == ComponentType.XZERO:
                    label = "X0"
            else:
                if isinstance(comp.type, ComponentType):
                    label = f"{comp.type.name}: {comp_name}"
                else:
                    label = f"{comp.type}: {comp_name}"

            dot.node(comp_name, label=label)

        # Render Bonds (Source --> Target)
        for i, bond in enumerate(bg.bonds):
            src_comp = bond.source.component.name
            tgt_comp = bond.target.component.name

            # Always direct edges from source to target
            dir_style = "forward"
            power_arrow = "halfopen"
            signal_arrow = "normal"

            if bond.target.causality == True and bond.type == ConnectionType.POWER_BOND:
                # Power arrow AND Causal stroke at target end
                arrowhead = f"tee{power_arrow}"
                arrowtail = "none"
            elif bond.source.causality == True and bond.type == ConnectionType.POWER_BOND:
                # Power arrow at target end, Causal stroke at source end
                arrowhead = power_arrow
                arrowtail = "tee"
                dir_style = "both"
            elif bond.type == ConnectionType.POWER_BOND:
                # Unassigned causality (only power flow arrow)
                arrowhead = power_arrow
                arrowtail = "none"
            else:
                # Signal bond (normal arrow)
                arrowhead = signal_arrow
                arrowtail = "none"
            if bond.type == ConnectionType.SIGNAL_BOND:
                label = f"s{i+1}"
            else:
                label = f"e{i+1}, f{i+1}"
            dot.edge(
                src_comp,
                tgt_comp,
                label=label,
                fontname="Helvetica-Oblique",
                fontsize="11",
                dir=dir_style,
                arrowhead=arrowhead,
                arrowtail=arrowtail,
                arrowsize="1.0",
                penwidth="1.5"
            )

        dot.render(filename=filename, format=format, cleanup=True, view=view)
        return dot