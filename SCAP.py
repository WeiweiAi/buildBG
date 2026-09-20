from enum import Enum, auto
from defineBG import  ComponentType,  Component, Bond, BondGraph, StorageType, ConnectionType, exportBG,importBG
from collections import deque

# System classifications returned after causality assignment completes.
class SystemType(Enum):
    """Classifies a graph's resulting ordinary or differential-algebraic system."""
    ODE = auto()
    DAE_DERIVATIVE = auto()
    DAE_ALGEBRAIC = auto()
    DAE_MIXED = auto()
class SCAPEngine:
    """Sequential Causality Assignment Procedure."""
    def __init__(self, graph: BondGraph):
        """Creates a causality solver that operates directly on `graph`."""
        self.graph = graph

    def _get_bond_effort_direction(self, bond: Bond, component: Component) -> str | None:
        """Determines if effort is flowing 'IN' to or 'OUT' of the given component via this bond."""
        if bond.source.causality is None or bond.target.causality is None:
            return None # Causality not yet assigned for this bond
            
        if bond.source.component == component:
            return "IN" if bond.source.causality  else "OUT"
        elif bond.target.component == component:
            return "IN" if bond.target.causality  else "OUT"
        return None 

    def _propagate_worklist(self, initial_components: list[Component]) -> None:
        """Propagates causality constraints through queued neighboring components."""
        # Avoid duplicate queue entries: a component waiting in the queue
        # does not need to be scheduled again.
        worklist = deque()
        queued: set[Component] = set()

        for comp in initial_components:
            if comp not in queued:
                worklist.append(comp)
                queued.add(comp)

        while worklist:
            comp = worklist.popleft()
            queued.discard(comp)
            self._propagate_component_causality(comp, worklist, queued)

    def _propagate_component_causality(self,comp: Component, worklist: deque[Component], queued: set[Component] ) -> bool:
        """Applies the local junction or transducer causality rules for one component."""
        
        assigned_bonds: list[Bond] = []
        unassigned_bonds: list[Bond] = []
        for bond in self.graph._get_bonds_for_component(comp):
            if bond.source.causality is None or bond.target.causality is None:
                unassigned_bonds.append(bond)
            else:
                assigned_bonds.append(bond)

        if not unassigned_bonds:
            return False

        changed_bonds: list[Bond] = []

        # 0-Junction & Switched 0-Junction (1 Effort IN constraint)
        if comp.type in (ComponentType.ZERO, ComponentType.XZERO):
            effort_in = sum(1 for b in assigned_bonds if self._get_bond_effort_direction(b, comp) == "IN")       
            if effort_in > 1:
                raise ValueError(f"Causality Conflict: 0-Junction '{comp.name}' has {effort_in} effort inputs (max 1).")
            if effort_in == 1:
                for b in unassigned_bonds:
                    b.assign_causality(b.source, False) if b.source.component == comp else b.assign_causality(b.source, True)
                    changed_bonds.append(b)
            elif len(unassigned_bonds) == 1 and effort_in == 0:
                b = unassigned_bonds[0]
                b.assign_causality(b.source, False) if b.source.component == comp else b.assign_causality(b.source, True)
                changed_bonds.append(b)

        # 1-Junction & Switched 1-Junction (1 Effort OUT constraint)
        elif comp.type in (ComponentType.ONE, ComponentType.XONE):
            effort_out = sum(1 for b in assigned_bonds if self._get_bond_effort_direction(b, comp) == "OUT")

            if effort_out > 1:
                raise ValueError(f"Causality Conflict: 1-Junction '{comp.name}' has {effort_out} effort outputs / flow inputs (max 1).")

            if effort_out == 1:
                for b in unassigned_bonds:
                    b.assign_causality(b.source, True) if b.source.component == comp else b.assign_causality(b.source, False)
                    changed_bonds.append(b)
            elif len(unassigned_bonds) == 1 and effort_out == 0:
                b = unassigned_bonds[0]
                b.assign_causality(b.source, True) if b.source.component == comp else b.assign_causality(b.source, False)
                changed_bonds.append(b)

        # Transformers (TF, MTF)
        elif comp.type in (ComponentType.TF, ComponentType.MTF) and len(unassigned_bonds) == 1 and len(assigned_bonds) == 1:
            assigned_dir = self._get_bond_effort_direction(assigned_bonds[0], comp)
            b = unassigned_bonds[0]
            if assigned_dir is None:
                raise ValueError(f"Cannot determine effort direction for assigned bond on component '{comp.name}'.")
            else:
                if assigned_dir == "IN":
                    b.assign_causality(b.source, False) if b.source.component == comp else b.assign_causality(b.source, True)
                else:
                    b.assign_causality(b.source, True) if b.source.component == comp else b.assign_causality(b.source, False)
            changed_bonds.append(b)

        # Gyrators (GY, MGY)
        elif comp.type in (ComponentType.GY, ComponentType.MGY) and len(unassigned_bonds) == 1 and len(assigned_bonds) == 1:
            assigned_dir = self._get_bond_effort_direction(assigned_bonds[0], comp)
            b = unassigned_bonds[0]
            if assigned_dir is None:
                raise ValueError(f"Cannot determine effort direction for assigned bond on component '{comp.name}'.")
            else:
                if assigned_dir == "IN":
                    b.assign_causality(b.source, True) if b.source.component == comp else b.assign_causality(b.source, False)
                else:
                    b.assign_causality(b.source, False) if b.source.component == comp else b.assign_causality(b.source, True)
            changed_bonds.append(b)

        for b in changed_bonds:
            other_comp = b.target.component if b.source.component == comp else b.source.component
            if other_comp not in queued:
                worklist.append(other_comp)
                queued.add(other_comp)

        return len(changed_bonds) > 0

    def _trace_algebraic_loop(self, start_comp: Component, initial_bond: Bond) -> list[Bond]:
        """Traces an unbranched path and returns it when it closes at the start."""
        # Keep a set for O(1) visited-bond membership tests, while retaining
        # a list for the returned path.
        visited_bonds: set[Bond] = {initial_bond}
        loop_path: list[Bond] = [initial_bond]
        curr = (
            initial_bond.target.component
            if initial_bond.source.component == start_comp
            else initial_bond.source.component
        )

        while curr and curr != start_comp:
            next_bond = next(
                (b for b in self.graph._get_bonds_for_component(curr) if b not in visited_bonds),
                None,
            )
            if next_bond is None:
                break                                        
                                                                                                            
            visited_bonds.add(next_bond)
            loop_path.append(next_bond)
            curr = (
                next_bond.target.component
                if next_bond.source.component == curr
                else next_bond.source.component
            )

        return loop_path if curr == start_comp else []

    def _classify_system(self) -> SystemType:
        """Derives the system type from derivative causality and algebraic loops."""
        has_derivative = len(self.derivative_causality_components) > 0
        has_algebraic = len(self.algebraic_loops) > 0

        if has_derivative and has_algebraic:
            return SystemType.DAE_MIXED
        elif has_derivative:
            return SystemType.DAE_DERIVATIVE
        elif has_algebraic:
            return SystemType.DAE_ALGEBRAIC
        return SystemType.ODE
    
    def run(self) -> SystemType:
        """Assigns all power-bond causalities and classifies the resulting system."""

        for comp in self.graph.components.values():
            for port in comp.ports.values():
                port.causality = port.fixed_causality   

        self.derivative_causality_components = []
        self.algebraic_loops = []

        # =====================================================================
        # STEP 1: Fixed Causality Type 1a (Independent & Modulated Sources)
        # =====================================================================
        step123_neighbors: list[Component] = []
        for comp in self.graph.components.values():
            if comp.type in (ComponentType.SE, ComponentType.MSE, ComponentType.SF, ComponentType.MSF):
                for port in comp.ports.values():
                    if port.bond and port.bond.type == ConnectionType.POWER_BOND: # active bond
                        is_effort_source = comp.type in (ComponentType.SE, ComponentType.MSE)
                        if is_effort_source:
                            target_causality = False # provides effort and receives flow
                        else:
                            target_causality = True # receives effort and provides flow
                        port.bond.assign_causality(port, target_causality)
                        neighbor = port.bond.get_other_component(comp)
                        step123_neighbors.append(neighbor)
        # =====================================================================
        # STEP 2: Fixed Causality Type 1b (Non-Invertible / Blocks / Switched)
        # =====================================================================
        for comp in self.graph.components.values():
            # Targets explicit non-invertibles, signal blocks, or locked switches
            if getattr(comp, "non_invertible", False):
                for port in comp.ports.values():
                    if port.bond and port.bond.type == ConnectionType.POWER_BOND:
                        if port.fixed_causality is not None:
                            port.bond.assign_causality(port, port.fixed_causality)
                        else:
                            pass                       
                        neighbor = port.bond.get_other_component(comp)
                        step123_neighbors.append(neighbor)

        # =====================================================================
        # STEP 3: Preferred Causality (Integral Causality for Storage Elements)
        # =====================================================================
        storage_types = (ComponentType.C, ComponentType.MC, ComponentType.I, ComponentType.MI, ComponentType.IC, ComponentType.MIC)
        for comp in self.graph.components.values():
            if comp.type in storage_types:
                for port in comp.ports.values():
                    if port.bond and port.bond.type == ConnectionType.POWER_BOND:
                        # Determine if this specific port acts as C or I
                        # For mixed IC/MIC, check port-level definitions; fallback to component level
                        pref_causality = None
                        if port.storage == StorageType.C_TYPE:
                            pref_causality = False
                        elif port.storage == StorageType.I_TYPE:
                            pref_causality = True
                        if pref_causality is None:
                            continue # Skip if no preferred causality can be determined for this port
                        if port.causality is None:
                            port.bond.assign_causality(port, pref_causality)
                            neighbor = port.bond.get_other_component(comp)
                            step123_neighbors.append(neighbor)
                        elif port.causality != pref_causality:
                            # Preferred causality was already lost.
                            # Preserve the existing assignment and make sure
                            # the opposite endpoint is assigned consistently.
                            port.bond.assign_causality(port, port.causality)
                            if comp not in self.derivative_causality_components:
                                self.derivative_causality_components.append(comp)

        self._propagate_worklist(step123_neighbors)
        # =====================================================================
        # STEP 4: Arbitrary / Free Causality & Algebraic Loop Inventory
        # =====================================================================
        for comp in self.graph.components.values():
            if comp.type in (ComponentType.R, ComponentType.MR):
                for port in comp.ports.values():
                    if port.bond and port.bond.type == ConnectionType.POWER_BOND:
                        # Assign arbitrary effort out
                        if port.causality is None:
                            if other_port := port.bond.get_other_port(port):
                                if other_port.causality is None:
                                    port.bond.assign_causality(port, True) # Effort at source, flow at target
                                else:
                                    port.bond.assign_causality(port, not other_port.causality)  
                        else:
                            port.bond.assign_causality(port, port.causality) # Effort at source, flow at target                      
                        # Trace if this choice formed an algebraic loop back to itself
                        loop_bonds = self._trace_algebraic_loop(comp, port.bond)
                        if loop_bonds:
                            self.algebraic_loops.append(loop_bonds)

                        neighbor = port.bond.get_other_component(comp)
                        self._propagate_worklist([neighbor])

        unassigned = [b for b in self.graph.bonds if b.type == ConnectionType.POWER_BOND and (b.source.causality is None or b.target.causality is None)]
        if unassigned:
            print(f"Unassigned bonds: {[b.name for b in unassigned]}")
            raise RuntimeError(f"SCAP failed: {len(unassigned)} bond(s) remained unassigned. Check ill-posed structures or disconnected loops.")
        for bond in self.graph.bonds:
            bond.validate_causality() # Ensure all bonds are valid after causality assignment
        self.system_type = self._classify_system()
        return self.system_type

if __name__ == "__main__":
    bg = importBG('mass_spring_damper.json')  # Load a BondGraph from a JSON file
    # Add components and bonds to the bond graph as needed
    scap_engine = SCAPEngine(bg)
    system_type = scap_engine.run()
    exportBG(bg, 'mass_spring_damper_causality.json')  # Export the BondGraph to a JSON file
    print(f"System type: {system_type.name}")  
    