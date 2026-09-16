from defineBG import JUNCTIONS, Component, Port, Bond, ComponentType, ConnectionType, BondGraph, Equation, importBG,exportBG
import sympy as sp
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
    def generate_network_equations(self) -> list[Equation]:
        """
        Generates network equations by solving the entire topological system simultaneously,
        eliminating all internal junction variables.
        """
        for bond in self.graph.bonds:
            try:
                bond.validate_causality()
            except ValueError as e:
                raise RuntimeError(f"Cannot generate equations: {e} due to non-deterministic causality assignment.") from e

        system_equations = []
        dependent_symbols = []
        junction_symbols = [] # Track all junction variables to eliminate them
        port_symbol_map = {}

        # 1. Bond Continuity Equations
        for bond in self.graph.bonds:
            if bond.type == ConnectionType.POWER_BOND:
                e_src, e_tgt = sp.Symbol(bond.source.effort.id), sp.Symbol(bond.target.effort.id)
                f_src, f_tgt = sp.Symbol(bond.source.flow.id), sp.Symbol(bond.target.flow.id)
                system_equations.append(sp.Eq(e_src, e_tgt))
                system_equations.append(sp.Eq(f_src, f_tgt))

        # 2. Junction Laws (and tracking junction symbols)
        for comp in self.graph.components.values():
            if comp.type in (ComponentType.ZERO, ComponentType.XZERO):
                flows_in = []
                flows_out = []
                base_effort = None
                
                for bond in comp.bonds:
                    if bond.type != ConnectionType.POWER_BOND: continue
                    port = self._get_port_for_component(bond, comp)
                    
                    # Track junction symbols so SymPy knows to eliminate them
                    junction_symbols.append(sp.Symbol(port.effort.id))
                    junction_symbols.append(sp.Symbol(port.flow.id))
                    
                    if base_effort is None:
                        base_effort = sp.Symbol(port.effort.id)
                    else:
                        system_equations.append(sp.Eq(base_effort, sp.Symbol(port.effort.id)))
                        
                    sym = sp.Symbol(port.flow.id)
                    if port is bond.target:
                        flows_in.append(sym)
                    else:
                        flows_out.append(sym)
                
                if base_effort is not None:
                    system_equations.append(sp.Eq(sp.Add(*flows_in), sp.Add(*flows_out)))

            elif comp.type in (ComponentType.ONE, ComponentType.XONE):
                efforts_in = []
                efforts_out = []
                base_flow = None
                
                for bond in comp.bonds:
                    if bond.type != ConnectionType.POWER_BOND: continue
                    port = self._get_port_for_component(bond, comp)
                    
                    # Track junction symbols so SymPy knows to eliminate them
                    junction_symbols.append(sp.Symbol(port.effort.id))
                    junction_symbols.append(sp.Symbol(port.flow.id))
                    
                    if base_flow is None:
                        base_flow = sp.Symbol(port.flow.id)
                    else:
                        system_equations.append(sp.Eq(base_flow, sp.Symbol(port.flow.id)))
                        
                    sym = sp.Symbol(port.effort.id)
                    if port is bond.target:
                        efforts_in.append(sym)
                    else:
                        efforts_out.append(sym)
                
                if base_flow is not None:
                    system_equations.append(sp.Eq(sp.Add(*efforts_in), sp.Add(*efforts_out)))

            else:
                # 3. Component Dependent Variables
                for port in comp.ports.values():
                    if port.bond is None: continue                    
                    if port.causality is True:
                        dep_sym = sp.Symbol(port.effort.id)
                        dependent_symbols.append(dep_sym)
                    elif port.causality is False:
                        dep_sym = sp.Symbol(port.flow.id)
                        dependent_symbols.append(dep_sym)

        # 4. Solve the system
        equations = []
        try:
            # Solve for BOTH dependent variables AND junction variables
            # This forces SymPy to express the answers ONLY in terms of independent component variables
            symbols_to_solve = dependent_symbols + junction_symbols
            solutions = sp.solve(system_equations, symbols_to_solve, dict=True)            
            if solutions:
                sol_dict = solutions[0]
                # Filter the results to only output the equations for the components, discarding junction equations
                for dep_sym in dependent_symbols:
                    if dep_sym in sol_dict:
                        solved_expr = sol_dict[dep_sym]
                        equations.append(
                            Equation(
                                dependent_symbol=f"{dep_sym}",
                                infix_rhs=f"{solved_expr}",
                                voi="",  # This could be set to the time variable if needed
                                expression=f"{dep_sym} = {solved_expr}",
                                description=f"Network equation"
                            )
                        )
        except Exception as e:
            print(f"Warning: Failed to solve network equations. Error: {e}")

        return equations

    def _process_ode_expression(self, expr: sp.Expr) -> Equation | None:
        ode_func = sp.Function('ode')
        found_odes = expr.find(ode_func)

        if not found_odes:
            return None 

        ode_term = list(found_odes)[0]

        # EXTRACT INNER ARGUMENTS: 
        # ode_term.args[0] is the state variable (e.g., x_C_Spring)
        # ode_term.args[1] is the time variable (e.g., t)
        state_variable = str(ode_term.args[0]) 
        voi_sym = str(ode_term.args[1])

        solutions = sp.solve(expr, ode_term)

        if not solutions:
            raise ValueError(f"Could not solve expression for {ode_term}")

        rhs_expr = solutions[0]
        infix_rhs_str = str(rhs_expr)

        formatted_expression = f"ode({state_variable}, {voi_sym}) = {infix_rhs_str}"

        return Equation(
            dependent_symbol=state_variable,  # Use the extracted inner variable here!
            infix_rhs=infix_rhs_str,
            voi=voi_sym,
            expression=formatted_expression,
            description=""
        )
    
    def _generate_constitutive_equations(self) -> list[Equation]:
        """
        Uses SymPy to solve component zero-form equations based on causality.
        """
        equations = []
        ode_func = sp.Function('ode') # Define this at the top

        subs_map_graph = {}
        if hasattr(self.graph, 'physical_constants') and self.graph.physical_constants is not None:
            for pq in self.graph.physical_constants:
                param_generic = sp.Symbol(pq.physical_quantity.symbol if pq.physical_quantity and pq.physical_quantity.symbol else pq.id)
                param_id = sp.Symbol(pq.id)
                subs_map_graph[param_generic] = param_id
        
        for comp in self.graph.components.values():
            if comp.type in JUNCTIONS:
                continue
                
            if not hasattr(comp, 'constitutive_equations') or not comp.constitutive_equations:
                continue
                
            dependent_symbols = []
            subs_map = {}
            if comp.parameters is not None:
                for pq in comp.parameters:
                    param_generic = sp.Symbol(pq.physical_quantity.symbol if pq.physical_quantity and pq.physical_quantity.symbol else pq.id)
                    param_id = sp.Symbol(pq.id)
                    subs_map[param_generic] = param_id
                    
            for i, (port_label, port) in enumerate(comp.ports.items()):
                e_generic, f_generic = sp.Symbol(f"e_{i}"), sp.Symbol(f"f_{i}")
                q_generic, p_generic = sp.Symbol(f"q_{i}"), sp.Symbol(f"p_{i}")
                s_generic = sp.Symbol(f"s_{i}")
                
                e_id, f_id = sp.Symbol(port.effort.id), sp.Symbol(port.flow.id)
                q_id, p_id = sp.Symbol(port.quantity.id), sp.Symbol(port.momentum.id)
                s_id = sp.Symbol(port.signal.id)
                
                subs_map.update({
                    e_generic: e_id, f_generic: f_id,
                    q_generic: q_id, p_generic: p_id, s_generic: s_id
                })
                
                # ONLY append algebraic variables (effort/flow) to dependent_symbols
                if port.causality is True:
                    dependent_symbols.append(f_id)
                elif port.causality is False:
                    dependent_symbols.append(e_id)

            # 3. Parse strings into SymPy objects (safely passing the ode_func)
            raw_exprs = [sp.parse_expr(eq, local_dict={'ode': ode_func}) for eq in comp.constitutive_equations]
            id_exprs = [expr.subs(subs_map | subs_map_graph) for expr in raw_exprs]
            
            # --- NEW LOGIC: Separate ODEs from Algebraic Constraints ---
            ode_exprs = []
            alg_exprs = []
            for expr in id_exprs:
                if expr.has(ode_func):
                    ode_exprs.append(expr)
                else:
                    alg_exprs.append(expr)

            # 4. Process the ODEs using your custom method
            for ode_expr in ode_exprs:
                eq = self._process_ode_expression(ode_expr)
                if eq:
                    eq.description = f"State equation for {comp.name}"
                    equations.append(eq)

            # 5. Process Algebraic equations using standard solve
            if alg_exprs and dependent_symbols:
                try:
                    solutions = sp.solve(alg_exprs, dependent_symbols, dict=True)
                    if solutions:
                        sol_dict = solutions[0]
                        for dep_sym, solved_expr in sol_dict.items():
                            equations.append(
                                Equation(
                                    dependent_symbol=f"{dep_sym}",
                                    infix_rhs=f"{solved_expr}",
                                    voi="",
                                    expression=f"{dep_sym} = {solved_expr}",
                                    description=f"{comp.type.name} constitutive equation for {comp.name}" if isinstance(comp.type, ComponentType) else f"Constitutive equation for {comp.name}"
                                )
                            )
                except Exception as e:
                    print(f"Warning: Could not solve algebraic equations for {comp.name}. Error: {e}")
                
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

                if eq.voi:
                    new_voi = sp.Symbol(eq.voi).subs(id_to_symbol)
                    new_dependent_symbol = sp.Symbol(eq.dependent_symbol).subs(id_to_symbol)
                else:
                    new_voi = ""
                    new_dependent_symbol = new_lhs
                translated_equations.append(
                    Equation(
                        dependent_symbol=f"{new_dependent_symbol}",
                        infix_rhs=f"{new_rhs}",
                        voi=f"{new_voi}",
                        expression=f"{new_lhs} = {new_rhs}",
                        description=eq.description
                    )
                )
            except Exception as e:
                print(f"Warning: Failed to translate equation '{eq.expression}'. Error: {e}")
                translated_equations.append(eq) # Return untranslated on failure

        self.graph.equations = translated_equations        
        return translated_equations    

if __name__ == "__main__":
    # Example usage
    bg = importBG("mass_spring_damper_causality.json")
    builder = EquationBuilder(bg)
    
    equations = builder.generate_network_equations()
    print("Network Equations:")
    for eq in equations:
        print("expression:", eq.expression, "\n","description:", eq.description)

    constitutive_eqs = builder._generate_constitutive_equations()
    print("\nConstitutive Equations:")
    for eq in constitutive_eqs:
        print("expression:", eq.expression, "\n","description:", eq.description)

    equations_with_symbols = builder.translate_to_symbols(equations + constitutive_eqs)
    print("\nEquations with Physical Symbols:")
    for eq in equations_with_symbols:
        print("dependent_symbol:", eq.dependent_symbol, "expression:", eq.expression, "\n","description:", eq.description)    

    exportBG(bg, "mass_spring_damper_equations.json")      