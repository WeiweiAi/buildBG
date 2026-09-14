from dataclasses import dataclass
from defineBG import JUNCTIONS,BGVariable, Component, Port, Bond, ComponentType, ConnectionType, BondGraph,importBG
@dataclass
class Equation:
    """Represents a single equation: y = f(x) with a description."""
    port_name: str = ""
    variable: BGVariable = BGVariable.EFFORT
    expression: str = ""  # Could be a string for now, or a sympy.Expr in a real solver
    description: str = ""
@dataclass(frozen=True)
class JunctionPropagation:
    determining_junction_port: Port # The port on the junction determining the effort (zero junction) or flow (one junction) values.
    determining_component_port: Port # The port of a component connected to determining_junction_port
    propagated_expression: str
    changed: bool
@dataclass(frozen=True)
class JunctionRule:
    propagated_variable: str
    conserved_variable: str
    determining_causality: bool

ZERO_RULE = JunctionRule(
    propagated_variable="effort",
    conserved_variable="flow",
    determining_causality=True,
)

ONE_RULE = JunctionRule(
    propagated_variable="flow",
    conserved_variable="effort",
    determining_causality=False,
) 

class EquationBuilder:
    """Mathematical model generation."""
    def __init__(self, graph: BondGraph):
        self.graph = graph

    def _get_port_for_component(
            self,
            bond: Bond,
            component: Component,
        ) -> Port:
            if bond.source.component is component:
                return bond.source

            if bond.target.component is component:
                return bond.target

            raise RuntimeError(
                f"Bond '{bond.name}' is not connected to "
                f"component '{component.name}'."
            )
    def _initialize_component_variables(self, variable_expr: dict) -> None:

         # ------------------------------------------------------------------
            # Initialize component-side known variables.
            #
            # True  -> receives effort -> flow is component-side known
            # False -> provides effort -> effort is component-side known
        # ------------------------------------------------------------------       

        for comp in self.graph.components.values():
            if comp.component_type in JUNCTIONS:
                continue
            for port in comp.ports.values():
                if port.bond is None:
                    continue

                if port.causality is True:
                    # Component receives effort and provides flow.
                    variable_expr[port]["flow"] = port.flow.id

                elif port.causality is False:
                    # Component receives flow and provides effort.
                    variable_expr[port]["effort"] = port.effort.id

                else:
                    raise RuntimeError(
                        f"Port '{port.name}' has no assigned causality."
                    )
    def _set_variable_expression(
            self,
            variable_expr: dict,
            port: Port,
            variable: str,
            expression: str,
        ) -> bool:
            """Set an expression and report whether it was newly established."""

            current = variable_expr[port][variable]

            if current is None:
                variable_expr[port][variable] = expression
                return True # Successfully set a new expression

            if current != expression:
                raise ValueError(
                    f"Conflicting {variable} expressions for "
                    f"port '{port.name}': "
                    f"'{current}' versus '{expression}'."
                )

            return False # Expression was already set and matches the new one
    
    def _propagate_junction(
            self,
            junction: Component,
            rule: JunctionRule,
            variable_expr: dict,
        ) -> JunctionPropagation | None:

            power_bonds = [bond
                for bond in junction.bonds
                if bond.connection_type == ConnectionType.POWER_BOND
            ]
            if len(power_bonds) < 2:
                return None

            junction_ports = [self._get_port_for_component(bond, junction)
                for bond in power_bonds
            ]
            source_ports = [ port 
                for port in junction_ports
                if port.causality is rule.determining_causality
            ]

            if len(source_ports) != 1:
                raise ValueError(
                    f"Junction '{junction.name}' must have exactly one "
                    f"causal source port with causality="
                    f"{rule.determining_causality}; found {len(source_ports)}."
                )

            source_port = source_ports[0]

            if source_port.bond is None:
                raise RuntimeError(
                    f"Junction port '{source_port.name}' has no bond."
                )

            source_other_port = source_port.bond.get_other_port(source_port)

            expression = variable_expr[source_other_port][rule.propagated_variable]

            if expression is None:
                return None

            changed = False

            # The junction-side variable is not needed in the final equations,
            # but keeping it in the propagation map is useful internally.
            if variable_expr[source_port][rule.propagated_variable] is None:
                variable_expr[source_port][rule.propagated_variable] = expression
                changed = True

            for port in junction_ports:
                if port is source_port:
                    continue
                
                bond = port.bond
                if bond is None:
                    continue
                
                other_port = bond.get_other_port(port)

                if self._set_variable_expression(
                    variable_expr,
                    other_port,
                    rule.propagated_variable,
                    expression,
                ):
                    changed = True

            return JunctionPropagation(
                determining_junction_port=source_port,
                determining_component_port=source_other_port,
                propagated_expression=expression,
                changed=changed,
            ) 
     
    def _generate_component_input_equations(self, variable_expr ) -> list[Equation]:

       # ------------------------------------------------------------------
         #  Generate network input equations for non-junction components.
         #
         # True  -> effort is network supplied
         # False -> flow is network supplied
        equations = []
        for comp in self.graph.components.values():            
            if comp.component_type in JUNCTIONS:
                continue                
            for port in comp.ports.values():                
                if port.bond is None:
                    continue
                
                if port.causality is True:
                    expr = variable_expr[port]["effort"]
                    if expr is None:
                       continue
                    equations.append(
                        Equation(f"{port.name}",
                                 BGVariable.EFFORT,
                            f"{port.effort.id} = {expr}",
                            f"Network effort equation for {port.name}",
                        )
                    )
                elif port.causality is False:
                    expr = variable_expr[port]["flow"]
                    if expr is None:
                        continue
                    equations.append(
                        Equation(
                            f"{port.name}",
                            BGVariable.FLOW,
                            f"{port.flow.id} = {expr}",
                            f"Network flow equation for {port.name}",
                        )
                    )
        return equations 
    def _junction_conservation(
            self,
            junction: Component,
            rule: JunctionRule,
            state: JunctionPropagation,
        ) -> Equation | None:
            """
            Generate the junction conservation equation in solved form.

            The dependent variable is determined by the source_port already
            identified during propagation.

            No junction-port symbols are used in the final equation.
            """

            terms: list[tuple[int, str, Port]] = []

            for bond in junction.bonds:
            
                if bond.connection_type != ConnectionType.POWER_BOND:
                    continue
                
                junction_port = self._get_port_for_component(
                    bond,
                    junction,
                )

                other_port = bond.get_other_port(junction_port)

                if junction_port is bond.target:
                    sign = +1
                elif junction_port is bond.source:
                    sign = -1
                else:
                    raise RuntimeError(
                        f"Bond '{bond.name}' is inconsistent with "
                        f"junction '{junction.name}'."
                    )

                expression = getattr(getattr( other_port, rule.conserved_variable), "id")

                terms.append(
                    (
                        sign,
                        expression,
                        other_port,
                    )
                )

            if len(terms) < 2:
                return None

            # ---------------------------------------------------------------
            # Reuse the source port identified during propagation.
            # Its opposite port is the dependent non-junction variable.
            # ---------------------------------------------------------------
            dependent_port = state.determining_component_port

            dependent_term = next(
                (
                    (sign, expression, port)
                    for sign, expression, port in terms
                    if port is dependent_port
                ),
                None,
            )

            if dependent_term is None:
                raise RuntimeError(
                    f"Dependent port '{dependent_port.name}' was not found "
                    f"in the conservation terms for junction "
                    f"'{junction.name}'."
                )

            dependent_sign, dependent_expression, _ = dependent_term

            # ---------------------------------------------------------------
            # Solve:
            #
            #     s_d*x_d + Σ s_i*x_i = 0
            #
            # for x_d:
            #
            #     x_d = -1/s_d * Σ s_i*x_i
            #
            # Since s_d = ±1:
            #
            #     x_d = Σ (-s_i*s_d)*x_i
            # ---------------------------------------------------------------
            rhs_terms: list[str] = []

            for sign, expression, port in terms:
            
                if port is dependent_port:
                    continue
                
                rhs_sign = -sign * dependent_sign

                if not rhs_terms:
                    prefix = "" if rhs_sign > 0 else "-"
                else:
                    prefix = " + " if rhs_sign > 0 else " - "

                rhs_terms.append(
                    f"{prefix}{expression}"
                )

            rhs = "".join(rhs_terms)

            if not rhs:
                rhs = "0"

            junction_type_name = junction.component_type.name if isinstance(junction.component_type, ComponentType) else str(junction.component_type)
            return Equation(f"{dependent_port.name}",
                            BGVariable[rule.conserved_variable.upper()],
                f"{dependent_expression} = {rhs}",
                f"{junction_type_name} "
                f"{rule.conserved_variable} conservation of {junction.name}",
            )
    def _generate_junction_equations(
            self,
            junction: Component,
            variable_expr: dict,
        ) -> list[Equation]:

            if junction.component_type in (
                ComponentType.ZERO,
                ComponentType.XZERO,
            ):
                rule = ZERO_RULE

            elif junction.component_type in (
                ComponentType.ONE,
                ComponentType.XONE,
            ):
                rule = ONE_RULE

            else:
                return []

            state = self._propagate_junction(
                junction,
                rule,
                variable_expr,
            )

            if state is None:
                return []

            equation = self._junction_conservation(
                junction,
                rule,
                state,
            )

            return [equation] if equation is not None else []
    
    
    def generate_network_equations(self) -> list[Equation]:
            """
            Generate the network equations after SCAP causality assignment.

            Workflow
            --------
            1. Initialize component-side known effort/flow variables.
            2. Propagate the common variable through 0/1 junctions.
            3. Generate input equations for non-junction components.
            4. Generate solved-form conservation equations for junctions.

            Causality convention
            --------------------
            port.causality is True:
                port receives effort
                → component provides flow
                → effort is supplied by the network

            port.causality is False:
                port provides effort
                → component receives flow
                → flow is supplied by the network

            Junction equations do not contain junction-port symbols.
            """
            # ------------------------------------------------------------------
            # 0. Require completed causality assignment
            # ------------------------------------------------------------------
            for bond in self.graph.bonds:
                try:
                    bond.validate_causality()
                except ValueError as e:
                    raise RuntimeError(f"Cannot generate equations: {e} due to non-deterministic causality assignment.") from e
            # ------------------------------------------------------------------
            # 1. Create the propagation map.
            # None means:
            #   "this variable has not been determined by network propagation."
            # ------------------------------------------------------------------
            variable_expr: dict[Port, dict[str, str | None]] = {
                port: {
                    "effort": None,
                    "flow": None,
                }
                for comp in self.graph.components.values()
                for port in comp.ports.values()
                if port.bond is not None
            }
           
            self._initialize_component_variables(variable_expr)

            #  Generate junction conservation equations.
            equations: list[Equation] = []
            for junction in self.graph.components.values():
                if junction.component_type not in JUNCTIONS:
                    continue            
                equation = self._generate_junction_equations(junction, variable_expr)                
                if len(equation) > 0 :
                    equations.extend(equation)
            # Generate input equations for non-junction components.
            equations=equations+self._generate_component_input_equations(variable_expr)
            return equations     

if __name__ == "__main__":
    # Example usage
    bg = importBG("mass_spring_damper_causality.json")
    builder = EquationBuilder(bg)
    equations = builder.generate_network_equations()
    for eq in equations:
        print("port_name:", eq.port_name, "variable:", eq.variable, "\n", "expression:", eq.expression, "\n","description:", eq.description)