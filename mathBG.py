from dataclasses import dataclass
from defineBG import JUNCTIONS,BGVariable, Component, Port, Bond, ComponentType, ConnectionType, BondGraph,importBG
import sympy as sp

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
            if comp.type in JUNCTIONS:
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
                if bond.type == ConnectionType.POWER_BOND
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
            if comp.type in JUNCTIONS:
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
            
                if bond.type != ConnectionType.POWER_BOND:
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

            junction_type_name = junction.type.name if isinstance(junction.type, ComponentType) else str(junction.type)
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

            if junction.type in (
                ComponentType.ZERO,
                ComponentType.XZERO,
            ):
                rule = ZERO_RULE

            elif junction.type in (
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
                if junction.type not in JUNCTIONS:
                    continue            
                equation = self._generate_junction_equations(junction, variable_expr)                
                if len(equation) > 0 :
                    equations.extend(equation)
            # Generate input equations for non-junction components.
            equations=equations+self._generate_component_input_equations(variable_expr)
            return equations 
    
    def _generate_constitutive_equations(self) -> list[Equation]:
        """
        Uses SymPy to solve component zero-form equations based on causality.
        """
        equations = []

        subs_map_graph = {}
        if hasattr(self.graph, 'physical_constants') and self.graph.physical_constants is not None:
            for pq in self.graph.physical_constants:
                param_generic = sp.Symbol(pq.physical_quantity.symbol if pq.physical_quantity and pq.physical_quantity.symbol else pq.id)
                param_id = sp.Symbol(pq.id)
                subs_map_graph[param_generic] = param_id
        
        for comp in self.graph.components.values():
            # Skip junctions, they are handled by the network equations
            if comp.type in JUNCTIONS:
                continue
                
            if not hasattr(comp, 'constitutive_equations') or not comp.constitutive_equations:
                continue
                
            dependent_symbols = []
            symbol_to_port = {}
            subs_map = {}
            if comp.parameters is not None:
                for pq in comp.parameters:
                    param_generic = sp.Symbol(pq.physical_quantity.symbol if pq.physical_quantity and pq.physical_quantity.symbol else pq.id)
                    param_id = sp.Symbol(pq.id)
                    subs_map[param_generic] = param_id
            # 1. Map generic variables (e_0, f_0) to unique port IDs
            for i, (port_label, port) in enumerate(comp.ports.items()):
                e_generic = sp.Symbol(f"e_{i}")
                f_generic = sp.Symbol(f"f_{i}")
                q_generic = sp.Symbol(f"q_{i}")
                p_generic = sp.Symbol(f"p_{i}")
                s_generic = sp.Symbol(f"s_{i}")
                
                e_id = sp.Symbol(port.effort.id)
                f_id = sp.Symbol(port.flow.id)
                q_id = sp.Symbol(port.quantity.id)
                p_id = sp.Symbol(port.momentum.id)
                s_id = sp.Symbol(port.signal.id)
                
                subs_map[e_generic] = e_id
                subs_map[f_generic] = f_id
                subs_map[q_generic] = q_id
                subs_map[p_generic] = p_id
                subs_map[s_generic] = s_id
                
                # 2. Determine what to solve for based on causality
                if port.causality is True:
                    # Port receives effort -> Component provides flow
                    dependent_symbols.append(f_id)
                    symbol_to_port[f_id] = (port.name, BGVariable.FLOW)
                elif port.causality is False:
                    # Port provides effort -> Component receives flow
                    dependent_symbols.append(e_id)
                    symbol_to_port[e_id] = (port.name, BGVariable.EFFORT)

            # 3. Parse strings into SymPy objects and substitute generic variables with IDs
            raw_exprs = [sp.parse_expr(eq) for eq in comp.constitutive_equations]
            id_exprs = [expr.subs(subs_map | subs_map_graph) for expr in raw_exprs]
            
            # 4. Solve the system of equations for the dependent variables
            try:
                # dict=True returns a list of dictionaries mapping dependent symbols to expressions
                solutions = sp.solve(id_exprs, dependent_symbols, dict=True)
                
                if solutions:
                    sol_dict = solutions[0]
                    for dep_sym, solved_expr in sol_dict.items():
                        port_name, var_type = symbol_to_port[dep_sym]
                        equations.append(
                            Equation(
                                port_name=port_name,
                                variable=var_type,
                                expression=f"{dep_sym} = {solved_expr}",
                                description=f"{comp.type.name} constitutive equation for {comp.name}" if isinstance(comp.type, ComponentType) else f"Constitutive equation for {comp.name}"
                            )
                        )
            except Exception as e:
                print(f"Warning: Could not solve constitutive equations for {comp.name}. Error: {e}")
                
        return equations

    def translate_to_symbols(self, equations: list[Equation]) -> list[Equation]:
        """
        Uses SymPy to safely substitute abstract IDs with physical domain symbols.
        """
        id_to_symbol = {}
        
        # 1. Build the global translation dictionary
        for comp in self.graph.components.values():
            for port in comp.ports.values():
                for var_type in ["effort", "flow", "quantity", "momentum", "signal"]:
                    bg_var = getattr(port, var_type, None)
                    if bg_var and hasattr(bg_var, 'id'):
                        # If a physical symbol was assigned via the DomainRefiner, use it
                        if hasattr(bg_var, 'physical_quantity') and bg_var.physical_quantity and bg_var.physical_quantity.symbol:
                            # Append the port name to ensure symbols are unique (e.g., u_comp_p1) replace comp.p1 dot with underscore
                            sym = f"{bg_var.physical_quantity.symbol}"
                            id_to_symbol[sp.Symbol(bg_var.id)] = sp.Symbol(sym)
                        else:
                            # Fallback to the abstract ID if no physics are assigned
                            id_to_symbol[sp.Symbol(bg_var.id)] = sp.Symbol(bg_var.id)
                            
            # Map component-specific parameters
            if hasattr(comp, 'parameters') and comp.parameters is not None:
                for pq in comp.parameters:
                    sym = pq.physical_quantity.symbol if pq.physical_quantity and pq.physical_quantity.symbol else pq.id
                    # Make parameter symbol unique to the component to avoid clashing
                    id_to_symbol[sp.Symbol(pq.id)] = sp.Symbol(f"{sym}")
                    
        # Map global parameters (like R, T, F)
        if hasattr(self.graph, 'physical_constants') and self.graph.physical_constants is not None:
            for pq in self.graph.physical_constants:
                sym = pq.physical_quantity.symbol if pq.physical_quantity and pq.physical_quantity.symbol else pq.id
                id_to_symbol[sp.Symbol(pq.id)] = sp.Symbol(sym)

        # 2. Translate the equations safely using SymPy
        translated_equations = []
        for eq in equations:
            try:
                # Split into LHS and RHS for substitution
                lhs_str, rhs_str = eq.expression.split("=")
                lhs_expr = sp.parse_expr(lhs_str.strip())
                rhs_expr = sp.parse_expr(rhs_str.strip())
                
                new_lhs = lhs_expr.subs(id_to_symbol)
                new_rhs = rhs_expr.subs(id_to_symbol)
                
                translated_equations.append(
                    Equation(
                        port_name=eq.port_name,
                        variable=eq.variable,
                        expression=f"{new_lhs} = {new_rhs}",
                        description=eq.description
                    )
                )
            except Exception as e:
                print(f"Warning: Failed to translate equation '{eq.expression}'. Error: {e}")
                translated_equations.append(eq) # Return untranslated on failure
                
        return translated_equations    

if __name__ == "__main__":
    # Example usage
    bg = importBG("mass_spring_damper_causality.json")
    builder = EquationBuilder(bg)
    
    equations = builder.generate_network_equations()
    print("Network Equations:")
    for eq in equations:
        print("port_name:", eq.port_name, "variable:", eq.variable, "\n", "expression:", eq.expression, "\n","description:", eq.description)

    constitutive_eqs = builder._generate_constitutive_equations()
    print("\nConstitutive Equations:")
    for eq in constitutive_eqs:
        print("port_name:", eq.port_name, "variable:", eq.variable, "\n", "expression:", eq.expression, "\n","description:", eq.description)

    equations_with_symbols = builder.translate_to_symbols(equations + constitutive_eqs)
    print("\nEquations with Physical Symbols:")
    for eq in equations_with_symbols:
        print("port_name:", eq.port_name, "variable:", eq.variable, "\n", "expression:", eq.expression, "\n","description:", eq.description)    