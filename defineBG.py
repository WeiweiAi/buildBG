from __future__ import annotations
from enum import Enum, auto
from dataclasses import dataclass, field
from collections import deque
from typing import  Any,Callable
import warnings
import networkx as nx
import graphviz

# Sequential Causality Assignment Procedure (SCAP)
# Karnopp, Dean. "Alternative bond graph causal patterns and equation formulations for dynamic systems." (1983): 58-63.
class ComponentType(Enum):
    # Borutzky, Wolfgang. Bond graph modelling of engineering systems. Vol. 103. New York: springer, 2011.
    C = auto()  # Capacitor / Spring, a storage node only contains C-type ports, the flow variable is integrated with respect to time
    I = auto()  # Inductor / Mass,  a storage node only contains I-type ports, the effort variable is integrated with respect to time
    IC = auto() # storage node contains both C-type and I-type ports
    MC = auto() # A modulated C-type storage
    MI = auto() # A modulated I-type storage
    MIC = auto() # A modulated storage node contains both C-type and I-type ports
    R = auto()  # Resistor / Damper,https://bg-rdf.org/ontologies/bondgraph-framework#Dissipator
    MR = auto() # A modulated R-type dissipator
    # Sources
    SE = auto() # Effort Source, the dependent port variable is an effort
    SF = auto() # Flow Source, the dependent port variable is a flow
    MSE = auto() # Modulated Effort Source, the dependent port variable is an effort
    MSF = auto() # Modulated Flow Source, the dependent port variable is a flow
    SS = auto() # Source of Signal (for control systems)
    # Junctions
    ZERO = auto() # 0-junction (Common Effort)
    ONE = auto()  # 1-junction (Common Flow)
    XZERO = auto() # 0-junction with a boolean switch
    XONE = auto()  # 1-junction with a boolean switch
    # Transducers
    TF = auto()   # Transformer
    GY = auto()   # Gyrator
    MTF = auto()  # Modulated Transformer
    MGY = auto()  # Modulated Gyrator
    # User-defined / Custom
    CUSTOM = auto() # User-defined component
    BLOCK = auto()  # All nodes that only have signal ports and represent mathematical relations between these signals.

class PowerVariable(Enum):
    EFFORT = auto() # e.g., voltage, force, pressure
    FLOW = auto()   # e.g., current, velocity, volumetric flow rate
    QUANTITY = auto() # e.g., charge, displacement, volume
    MOMENTUM = auto() # e.g., momentum, angular momentum
    POWER = auto() # e.g., power, energy rate
    ENERGY = auto() # e.g., energy, work
    SIGNAL = auto() # A signal represents one arbitrary variable of time that may also be an effort or a flow, but not necessarily

class ConnectionType(Enum):
    POWER_BOND = auto() # carries both effort and flow (bidirectional)
    SIGNAL_BOND = auto() # carries only signal (unidirectional)

class PortType(Enum):
    POWER_PORT = auto() # Port for power exchange (effort and flow)
    SIGNAL_PORT = auto() # Port for signal interface (control signals)
    C_TYPE_PORT = auto() # Port for C-type storage (integrates flow to quantity)
    I_TYPE_PORT = auto() # Port for I-type storage (integrates effort to momentum)

class Causality(Enum):
    EFFORT_AT_SOURCE = auto() # Causal stroke at the source port, 
                              # i.e., the source port receives the effort (input) from the bond and the target port provides the effort (output) to the bond.
                              # i.e., the source port provides the flow (output) to the bond and the target port receives the flow (input) from the bond.
    EFFORT_AT_TARGET = auto() # Causal stroke at the target port, 
                              # i.e., the target port receives the effort (input) from the bond and the source port provides the effort (output) to the bond.
                              # i.e., the target port provides the flow (output) to the bond and the source port receives the flow (input) from the bond.
    UNASSIGNED = auto()

class ConstitutiveRelationship(Enum):
    PHI_C=auto() # q - PHI_C(e) = 0
    PHI_I=auto() # p - PHI_I(f) = 0
    PHI_IC=auto() # q - PHI_IC(e) = 0, p - PHI_IC(f) = 0
    PHI_R=auto() # e - PHI_R(f) = 0
    PHI_MC=auto() # q - PHI_MC(e) = 0
    PHI_MI=auto() # p - PHI_MI(f) = 0
    PHI_MIC=auto() # q - PHI_MIC(e) = 0, p - PHI_MIC(f) = 0
    PHI_MR=auto() # e - PHI_MR(f) = 0
    PHI_TF=auto() # e1 - PHI_TF(e2) = 0, f2 - PHI_TF(f1) = 0
    PHI_GY=auto() # e1 - PHI_GY(f2) = 0, e2 - PHI_GY(f1) = 0
    PHI_MTF=auto() # e1 - PHI_MTF(e2) = 0, f2 - PHI_MTF(f1) = 0
    PHI_MGY=auto() # e1 - PHI_MGY(f2) = 0, e2 - PHI_MGY(f1) = 0
    PHI_SE=auto() # e - PHI_SE(t) = 0
    PHI_SF=auto() # f - PHI_SF(t) = 0
    PHI_MSE=auto() # e - PHI_MSE(t) = 0
    PHI_MSF=auto() # f - PHI_MSF(t) = 0
    PHI_USER=auto() # User-defined constitutive relationship

class Domain(Enum):
    ABSTRACT = auto()  # Default state: uses e, f, p, q
    ELECTRICAL = auto()
    MECHANICAL_TRANSLATIONAL = auto()
    MECHANICAL_ROTATIONAL = auto()
    HYDRAULIC = auto()
    CHEMICAL = auto()
    THERMAL = auto()
    CUSTOM = auto()

class SystemType(Enum):
    ODE = auto()
    DAE_DERIVATIVE = auto()
    DAE_ALGEBRAIC = auto()
    DAE_MIXED = auto()

@dataclass
class PhysicalQuantity:
    """Metadata for a domain-specific power variable."""
    description: str
    symbol: str
    units: str

DOMAIN_VARIABLES: dict[Domain, dict[PowerVariable, PhysicalQuantity]] = {
    Domain.ABSTRACT: {
        PowerVariable.EFFORT: PhysicalQuantity("generalized effort", "e", "effort_units"),
        PowerVariable.FLOW: PhysicalQuantity("generalized flow", "f", "flow_units"),
        PowerVariable.QUANTITY: PhysicalQuantity("generalized extensive quantity", "q", "extensive_quantity_units"),
        PowerVariable.MOMENTUM: PhysicalQuantity("generalized momentum", "p", "momentum_units"),
    },
    Domain.ELECTRICAL: {
        PowerVariable.EFFORT: PhysicalQuantity("voltage", "u", "volt"),
        PowerVariable.FLOW: PhysicalQuantity("current", "i", "fA"),
        PowerVariable.QUANTITY: PhysicalQuantity("charge", "q", "fC"),
        PowerVariable.MOMENTUM: PhysicalQuantity("magnetic flux linkage", "p", "volt_s"),
    },
    Domain.MECHANICAL_TRANSLATIONAL: {
        PowerVariable.EFFORT: PhysicalQuantity("force", "F", "J_per_um"),
        PowerVariable.FLOW: PhysicalQuantity("velocity", "v", "um_per_s"),
        PowerVariable.QUANTITY: PhysicalQuantity("displacement", "x", "um"),
        PowerVariable.MOMENTUM: PhysicalQuantity("momentum", "p", "J_s_per_um"),
    },
    Domain.MECHANICAL_ROTATIONAL: {
        PowerVariable.EFFORT: PhysicalQuantity("torque", "T", "J_per_rad"),
        PowerVariable.FLOW: PhysicalQuantity("angular velocity", "w", "rad_per_s"),
        PowerVariable.QUANTITY: PhysicalQuantity("angular displacement", "theta", "rad"),
        PowerVariable.MOMENTUM: PhysicalQuantity("angular momentum", "p", "J_s_per_rad"),
    },
    Domain.HYDRAULIC: {
        PowerVariable.EFFORT: PhysicalQuantity("pressure", "P", "mmHg"),
        PowerVariable.FLOW: PhysicalQuantity("volume flow", "Q", "mL_per_s"),
        PowerVariable.QUANTITY: PhysicalQuantity("volume", "V", "mL"),
        PowerVariable.MOMENTUM: PhysicalQuantity("momentum of a flow tube", "p", "mmHg_mL2_per_s3"),
    },
    Domain.CHEMICAL: {
        PowerVariable.EFFORT: PhysicalQuantity("chemical potential", "mu", "J_per_mol"),
        PowerVariable.FLOW: PhysicalQuantity("molar flow", "v", "fmol_per_s"),
        PowerVariable.QUANTITY: PhysicalQuantity("molar amount", "q", "fmol")
    },
}

@dataclass
class StateVariable:
    """Represents a time-integrated energy state of a component (q or p)."""
    variable_type: PowerVariable
    component: Component = field(repr=False) # Prevents Infinite Recursion Crashing

    @property
    def symbol(self) -> str:
        """Dynamically resolves the physical symbol based on the domain."""
        domain_dict = DOMAIN_VARIABLES.get(self.component.domain)
        if domain_dict and self.variable_type in domain_dict:
            base_symbol = domain_dict[self.variable_type].symbol
            return f"{base_symbol}_{self.component.name}"
        
        # Fallback if domain is missing
        fallback = "q" if self.variable_type == PowerVariable.QUANTITY else "p"
        return f"{fallback}_{self.component.name}"
    @property
    def derivative_symbol(self) -> str:
        """Returns the time derivative of the state variable (x_dot)."""
        return f"d({self.symbol})/dt"

@dataclass
class ConstitutiveEquation:
    """Represents a single implicit relation: Phi(e, f, x, x_dot) = 0"""
    expression: Any  # Could be a string for now, or a sympy.Expr in a real solver
    description: str = ""

@dataclass
class Port:
    label: str
    component: Component = field(repr=False) # Prevents Infinite Recursion Crashing
    port_type: PortType = PortType.POWER_PORT
    fixed_causality: Causality | None = None # Either EFFORT_AT_SOURCE (output effort) or EFFORT_AT_TARGET (output flow) 
    domain: Domain = Domain.ABSTRACT
    bond: Bond | None = field(default=None, repr=False, init=False)

    @property
    def effective_domain(self) -> Domain:
        return self.domain if self.domain != Domain.ABSTRACT else self.component.domain

    @property
    def name(self) -> str:
        return f"{self.component.name}.{self.label}"
    
    def _get_symbol(self, variable_type: PowerVariable, default_prefix: str) -> str:
        domain_dict = DOMAIN_VARIABLES.get(self.effective_domain)
        if domain_dict and variable_type in domain_dict:
            return f"{domain_dict[variable_type].symbol}_{self.name}"
        return f"{default_prefix}_{self.name}"

    @property
    def effort(self) -> str:
        return self._get_symbol(PowerVariable.EFFORT, "e")

    @property
    def flow(self) -> str:
        return self._get_symbol(PowerVariable.FLOW, "f")

    @property
    def signal(self) -> str:
        return self._get_symbol(PowerVariable.SIGNAL, "s")

@dataclass
class Bond:
    source: Port
    target: Port
    connection_type: ConnectionType = ConnectionType.POWER_BOND
    causality: Causality = Causality.UNASSIGNED

    def __post_init__(self) -> None:
        # 1. Enforce the "exactly one bond per port" rule
        if self.source.bond is not None:
            raise ValueError(f"Source port {self.source.name} is already connected to a bond.")
        if self.target.bond is not None:
            raise ValueError(f"Target port {self.target.name} is already connected to a bond.")        
        # 2. If the checks pass, establish the bidirectional link
        self.source.bond = self
        self.target.bond = self

    @property
    def name(self) -> str:
        return f"{self.source.name}--{self.target.name}"

    @property
    def effort(self) -> str:
        if self.causality == Causality.EFFORT_AT_SOURCE: # Causal stroke at the source port
            return self.target.effort
        elif self.causality == Causality.EFFORT_AT_TARGET: # Causal stroke at the target port
            return self.source.effort
        else:
            if self.connection_type == ConnectionType.SIGNAL_BOND:
                # For signal bonds, we can default to source effort if causality is unassigned
                return self.source.effort
            else:
                raise ValueError("Causality is not assigned.")
        
    @property
    def flow(self) -> str:
        if self.causality == Causality.EFFORT_AT_SOURCE:
            return self.source.flow
        elif self.causality == Causality.EFFORT_AT_TARGET:
            return self.target.flow
        else:
            if self.connection_type == ConnectionType.SIGNAL_BOND:
                # For signal bonds, we can default to source flow if causality is unassigned
                return self.source.flow
            else:
                raise ValueError("Causality is not assigned.")
    @property
    def signal(self) -> str:
        if self.connection_type != ConnectionType.SIGNAL_BOND:
            raise ValueError("Signal property is only valid for signal bonds.")
        return self.source.signal

@dataclass(eq=False)
class Component:
    name: str
    component_type: ComponentType = ComponentType.CUSTOM
    initial_port_count: int = field(default=0, repr=False)
    domain: Domain = Domain.ABSTRACT
    non_invertible: bool = False # If True, the component has any constitutive relationship that cannot be algebraically inverted to solve for either effort or flow.
    ports: dict[str, Port] = field(default_factory=dict, repr=False, init=False)  
    # consider these later  
    parameters: dict[str, Any]  = field(default_factory=dict, repr=False, init=False)    
    states: dict[PowerVariable, StateVariable]  = field(default_factory=dict, repr=False, init=False)
    equations: list[ConstitutiveEquation] = field(default_factory=list, repr=False, init=False)
    # Optional callback to override the default linear equations
    equation_generator: Callable[[Component], list[ConstitutiveEquation]] | None = None

    def __hash__(self) -> int:
        # Hash only the unique string name, ignoring mutable fields like dicts
        return hash(self.name)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Component):
            return NotImplemented
        return self.name == other.name
    
    @property
    def port_count(self) -> int:
        """Dynamically always returns the true number of ports."""
        return len(self.ports)
    
    def __post_init__(self) -> None:
        # --- Flexible Auto-Port Generation ---
        # 1. Determine how many ports to generate        
        if self.initial_port_count==0: # if not specified
            # Fallback to standard defaults if not specified
            if self.component_type in (ComponentType.TF, ComponentType.GY, ComponentType.MC, ComponentType.MI, ComponentType.MR, ComponentType.MSE, ComponentType.MSF):
                self.initial_port_count = 2
            elif self.component_type in (ComponentType.MTF, ComponentType.MGY, ComponentType.MIC):
                self.initial_port_count = 3
            else:
                self.initial_port_count = 1 # by default, most components have at least one port
                
        # 2. Generate the assigned number of ports
        for i in range(1, self.initial_port_count + 1):
            label = f"p{i}"
            port_type = self._get_port_type(i)
            self.ports[label] = Port(label=label, component=self, port_type=port_type)

    def _get_port_type(self, port_index: int) -> PortType:
        """Determines the default port type based on the component type and port index."""
        if self.component_type in (ComponentType.C, ComponentType.MC):
            return PortType.C_TYPE_PORT if port_index ==1 else PortType.SIGNAL_PORT # 1st port is C-type, others are signal ports
        elif self.component_type in (ComponentType.I, ComponentType.MI):
            return PortType.I_TYPE_PORT if port_index ==1 else PortType.SIGNAL_PORT # 1st port is I-type, others are signal ports
        elif self.component_type in (ComponentType.IC, ComponentType.MIC):
            # For IC or MIC, we can alternate or assign based on index
            # 1st port is C-type, 2nd port is I-type, others are signal ports
            if port_index == 1:
                return PortType.C_TYPE_PORT
            elif port_index == 2:
                return PortType.I_TYPE_PORT
            else:
                return PortType.SIGNAL_PORT
        elif self.component_type in (ComponentType.R, ComponentType.MR):
            return PortType.POWER_PORT if port_index == 1 else PortType.SIGNAL_PORT # 1st port is power, others are signal ports
        elif self.component_type in (ComponentType.SE, ComponentType.SF, ComponentType.MSE, ComponentType.MSF):
            # 1 and 2nd ports are power ports, others are signal ports
            return PortType.POWER_PORT if port_index == 1 else PortType.SIGNAL_PORT
        elif self.component_type in (ComponentType.TF, ComponentType.GY, ComponentType.MTF, ComponentType.MGY):
            return PortType.POWER_PORT if port_index in (1, 2) else PortType.SIGNAL_PORT # 1st and 2nd ports are power, others are signal ports
        elif self.component_type == ComponentType.BLOCK:
            return PortType.SIGNAL_PORT
        else:
            return PortType.POWER_PORT  # Default to power port for other types

    def get_free_port(self) -> list[Port]:
        """
        Returns the available ports, or spawns a new one for junctions.
        Enforces port count limits for standard n-port elements.
        """
        # Junctions can dynamically spawn infinite ports
        if self.component_type in (ComponentType.ZERO, ComponentType.ONE):
            new_label = f"p{self.port_count + 1}"
            new_port = Port(label=new_label, component=self, port_type=PortType.POWER_PORT)
            self.ports[new_label] = new_port
            return [new_port]   
        # For all other elements, find an unconnected port
        return [port for port in self.ports.values() if port.bond is None]             
            
    def generate_atomic_equations(self) -> None: # Check this later
        self.equations.clear()
        
        if self.equation_generator is not None:
            self.equations.extend(self.equation_generator(self))
            return
            
        ports = list(self.ports.values())
        if not ports:
            return 
            
        # --- 1-Port Elements ---
        if self.component_type in (ComponentType.R, ComponentType.C, ComponentType.I, ComponentType.SE, ComponentType.SF):
            if len(ports) != 1:
                raise ValueError(f"{self.component_type.name} requires exactly 1 port.")
            p0 = ports[0]
            
            if self.component_type == ComponentType.R:
                self.equations.append(ConstitutiveEquation(f"{p0.effort} - ({self.parameters['R']}) * {p0.flow}", "Linear Dissipation"))
            elif self.component_type == ComponentType.C:
                state = self.states[PowerVariable.QUANTITY]
                self.equations.append(ConstitutiveEquation(f"{state.derivative_symbol} - {p0.flow}", "State Derivative"))
                self.equations.append(ConstitutiveEquation(f"{p0.effort} - (1/({self.parameters['C']})) * {state.symbol}", "Linear Energy Storage"))
            elif self.component_type == ComponentType.I:
                state = self.states[PowerVariable.MOMENTUM]
                self.equations.append(ConstitutiveEquation(f"{state.derivative_symbol} - {p0.effort}", "State Derivative"))
                self.equations.append(ConstitutiveEquation(f"{p0.flow} - (1/({self.parameters['I']})) * {state.symbol}", "Linear Energy Storage"))
            elif self.component_type == ComponentType.SE:
                self.equations.append(ConstitutiveEquation(f"{p0.effort} - SE_{self.name}(t)", "Effort Source"))
            elif self.component_type == ComponentType.SF:
                self.equations.append(ConstitutiveEquation(f"{p0.flow} - SF_{self.name}(t)", "Flow Source"))

        # --- 2-Port Transducers ---
        elif self.component_type in (ComponentType.TF, ComponentType.GY):
            if len(ports) < 2:
                raise ValueError(f"{self.component_type.name} requires at least 2 ports.")
            p0, p1 = ports[0], ports[1]
            
            if self.component_type == ComponentType.TF:
                m = self.parameters["m"]
                self.equations.append(ConstitutiveEquation(f"{p0.effort} - ({m}) * {p1.effort}", "TF Effort Relation"))
                self.equations.append(ConstitutiveEquation(f"({m}) * {p0.flow} - {p1.flow}", "TF Flow Relation"))   
            elif self.component_type == ComponentType.GY:
                r = self.parameters["r"]
                self.equations.append(ConstitutiveEquation(f"{p0.effort} - ({r}) * {p1.flow}", "GY Relation 1"))
                self.equations.append(ConstitutiveEquation(f"{p1.effort} - ({r}) * {p0.flow}", "GY Relation 2"))

        # --- Junctions (Rigorous Sign Convention) ---
        if self.component_type == ComponentType.ZERO:
            # 0-Junction: Common Effort (e1 = e2 = e3...)
            efforts = [p.effort for p in ports]
            for i in range(1, len(efforts)):
                self.equations.append(ConstitutiveEquation(f"{efforts[0]} - {efforts[i]}", "Common Effort"))               
            # Conservation of Flow: Sum(f_in) - Sum(f_out) = 0
            flow_terms = []
            for p in ports:
                if p.bond is None:
                    continue
                # If this port is the target of the bond, power flows INTO the junction (+)
                if p.bond.target == p:
                    flow_terms.append(f"+ {p.flow}")
                # If this port is the source of the bond, power flows OUT of the junction (-)
                else:
                    flow_terms.append(f"- {p.flow}")                   
            if flow_terms:
                sum_flows = " ".join(flow_terms).lstrip("+ ") # Clean up leading plus sign
                self.equations.append(ConstitutiveEquation(sum_flows, "Conservation of Flow"))

        if self.component_type == ComponentType.ONE:
            # 1-Junction: Common Flow (f1 = f2 = f3...)
            flows = [p.flow for p in ports]
            for i in range(1, len(flows)):
                self.equations.append(ConstitutiveEquation(f"{flows[0]} - {flows[i]}", "Common Flow"))
                
            # Conservation of Effort: Sum(e_in) - Sum(e_out) = 0
            effort_terms = []
            for p in ports:
                if p.bond is None:
                    continue
                # If this port is the target of the bond, power flows INTO the junction (+)
                if p.bond.target == p:
                    effort_terms.append(f"+ {p.effort}")
                # If this port is the source of the bond, power flows OUT of the junction (-)
                else:
                    effort_terms.append(f"- {p.effort}")
                    
            if effort_terms:
                sum_efforts = " ".join(effort_terms).lstrip("+ ") # Clean up leading plus sign
                self.equations.append(ConstitutiveEquation(sum_efforts, "Conservation of Effort"))


class BondGraph:
    def __init__(self, name: str = "bond_graph") -> None:
        self.name = name
        self.components: dict[str, Component] = {}
        self.bonds: list[Bond] = []

        # Extended Diagnostic State
        self.derivative_causality_components: list[Component] = []
        self.algebraic_loops: list[list[Bond]] = []
        self.system_type: SystemType = SystemType.ODE

    def to_networkx(self) -> nx.DiGraph:
        """Maps the completed bond graph structure to a NetworkX DiGraph."""
        G = nx.DiGraph(name=self.name)
        
        # Add components and ports
        for comp in self.components.values():
            G.add_node(comp.name, object=comp, kind="component")
            for port in comp.ports.values():
                G.add_node(port.name, object=port, kind="port")
                G.add_edge(comp.name, port.name, relationship="has_port")

        # Add bonds
        for bond in self.bonds:
            G.add_edge(
                bond.source.name, 
                bond.target.name, 
                object=bond, 
                kind="bond",
                causality=bond.causality.name
            )
            
        return G

    def add_component(self, component: Component | str, **kwargs) -> Component:
        """Adds a component to internal tracking."""
        if isinstance(component, str):
            comp_obj = Component(name=component, **kwargs)
        else:
            comp_obj = component

        if comp_obj.name in self.components:
            warnings.warn(f"Component '{comp_obj.name}' already exists.")
            return self.components[comp_obj.name]

        self.components[comp_obj.name] = comp_obj
        return comp_obj

    def add_port(self, comp_arg: Component | str, port_label: str, **kwargs) -> Port | None:
        """Creates a new port on a component, enforcing element port limits."""
        comp_name = comp_arg.name if isinstance(comp_arg, Component) else comp_arg
        comp = self.components.get(comp_name)
        
        if not comp:
            warnings.warn(f"Component '{comp_name}' not found in bond graph.")
            return None

        if port_label in comp.ports:
            warnings.warn(f"Port '{port_label}' already exists on '{comp.name}'. Returning existing port.")
            return comp.ports[port_label]

        if comp.component_type not in (ComponentType.ZERO, ComponentType.ONE):
            if len(comp.ports) >= comp.initial_port_count:
                warnings.warn(
                    f"Topology Error: Cannot add port '{port_label}' to '{comp.name}'. "
                    f"Component type {comp.component_type.name} is strictly limited to {comp.initial_port_count} port(s)."
                )
                return None

        new_port = Port(label=port_label, component=comp, **kwargs)
        comp.ports[port_label] = new_port
        return new_port

    def _auto_get_port(self, arg: Port | Component | str) -> list[Port]:
        """Resolves Port, Component, or String inputs into a valid Port object before adding a bond."""
        if isinstance(arg, Port):
            return [arg]

        if isinstance(arg, Component):
            return arg.get_free_port()

        if isinstance(arg, str):
            if "." in arg:
                comp_name, port_label = arg.split(".", 1)
                addPort=self.add_port(comp_name, port_label)
                if not addPort:
                    return []
                return [addPort]
            else:
                if arg in self.components:
                    return self._auto_get_port(self.components[arg])
                else:
                    warnings.warn(f"Component '{arg}' not found.")
                    return []

        warnings.warn(f"Cannot resolve bond endpoint from type: {type(arg)}")
        return []

    def add_bond(self, source: Component | Port | str, target: Component | Port | str, **kwargs) -> Bond | None:
        """Creates a bond between two endpoints."""
        src_port = self._auto_get_port(source)
        tgt_port = self._auto_get_port(target)

        if len(src_port) != 1 or len(tgt_port) != 1:
            warnings.warn("Each endpoint must resolve to exactly one port.")
            return None
        else:
            src_port = src_port[0]
            tgt_port = tgt_port[0]
            # check if source and target from the same component
            if src_port.component.name == tgt_port.component.name:
                warnings.warn("Cannot create a bond between ports of the same component.")
                return None
            bond = Bond(source=src_port, target=tgt_port, **kwargs)
            self.bonds.append(bond)
        return bond

    def _get_bonds_for_component(self, comp_arg: Component | str) -> list[Bond]:
        """Returns all bonds connected to a given component."""
        comp_name = comp_arg.name if isinstance(comp_arg, Component) else comp_arg
        comp = self.components.get(comp_name)
        if not comp:
            warnings.warn(f"Component '{comp_name}' not found in bond graph.")
            return []
        else:
            return [port.bond for port in comp.ports.values() if port.bond is not None]
    
    def _resolve_to_port(self, arg: Port | str) -> Port | None:
        """Resolves a Port, or String input into a valid Port object before deleting a bond."""
        if isinstance(arg, Port):
            return arg
        
        if isinstance(arg, str):
            if "." in arg:
                comp_name, port_label = arg.split(".", 1)
                comp = self.components.get(comp_name)
                if comp and port_label in comp.ports:
                    return comp.ports[port_label]
                else:
                    warnings.warn(f"Port '{port_label}' not found on component '{comp_name}'.")
                    return None
                
        warnings.warn(f"Cannot resolve to port from type: {type(arg)}")
        return None
    
    def delete_bond(self, arg1: Bond | Port | Component | str, arg2: Port | Component | str | None = None) -> int:
        """Deletes bonds based on flexible inputs."""
        bonds_to_delete: list[Bond] = []
        
        if arg2 is None:
            if not isinstance(arg1, Bond):
                warnings.warn("If only one argument is provided, it must be a Bond object.")
                return 0
            else:
                bonds_to_delete.append(arg1)
        else:
            if (isinstance(arg1, Port) or (isinstance(arg1, str) and "." in arg1)) and \
               (isinstance(arg2, Port) or (isinstance(arg2, str) and "." in arg2)):
                p1 = self._resolve_to_port(arg1)
                p2 = self._resolve_to_port(arg2)
                if not p1 or not p2:
                    warnings.warn("One or both ports could not be resolved. No bonds deleted.")
                    return 0
                if p1.bond and p2.bond and p1.bond == p2.bond:
                    bonds_to_delete.append(p1.bond)
                else:
                    warnings.warn("Ports are not connected by the same bond. No bonds deleted.")
                    return 0

            elif isinstance(arg1, (Component, str)) and isinstance(arg2, (Component, str)):
                bonds_c1 = self._get_bonds_for_component(arg1)
                bonds_c2 = self._get_bonds_for_component(arg2)
                bonds_to_delete = [bond for bond in bonds_c1 if bond in bonds_c2]
            else:
                warnings.warn("Invalid argument types for bond deletion. No bonds deleted.")
                return 0

        deleted_count = 0
        for bond in list(bonds_to_delete):
            if bond in self.bonds:
                bond.source.bond = None
                bond.target.bond = None
                self.bonds.remove(bond)
                deleted_count += 1

        return deleted_count

    def delete_component(self, comp_arg: Component | str) -> None:
        """Removes a component, its ports, and attached bonds from the graph."""
        comp_name = comp_arg.name if isinstance(comp_arg, Component) else comp_arg
        comp = self.components.get(comp_name)
        if not comp:
            warnings.warn(f"Component '{comp_name}' not found in bond graph.")
            return

        for bond in self._get_bonds_for_component(comp):
            self.delete_bond(bond)

        del self.components[comp.name]

    

    def _get_bond_effort_direction(self, bond: Bond, component: Component) -> str | None:
        """Determines if effort is flowing 'IN' to or 'OUT' of the given component via this bond."""
        if bond.causality == Causality.UNASSIGNED:
            return None
            
        if bond.source.component == component:
            return "OUT" if bond.causality == Causality.EFFORT_AT_TARGET else "IN"
        elif bond.target.component == component:
            return "IN" if bond.causality == Causality.EFFORT_AT_TARGET else "OUT"
        return None

    def assign_causality(self) -> SystemType:
        """Executes the Generalized Extended SCAP framework."""
        for bond in self.bonds:
            bond.causality = Causality.UNASSIGNED
            
        self.derivative_causality_components = []
        self.algebraic_loops = []

        # =====================================================================
        # STEP 1: Fixed Causality Type 1a (Independent & Modulated Sources)
        # =====================================================================
        step1_neighbors: list[Component] = []
        for comp in self.components.values():
            if comp.component_type in (ComponentType.SE, ComponentType.MSE, ComponentType.SF, ComponentType.MSF):
                for port in comp.ports.values():
                    if port.bond:
                        is_effort_source = comp.component_type in (ComponentType.SE, ComponentType.MSE)
                        target_causality = (
                            Causality.EFFORT_AT_TARGET if port.bond.source == port else Causality.EFFORT_AT_SOURCE
                        ) if is_effort_source else (
                            Causality.EFFORT_AT_SOURCE if port.bond.source == port else Causality.EFFORT_AT_TARGET
                        )

                        if port.bond.causality != Causality.UNASSIGNED and port.bond.causality != target_causality:
                            raise ValueError(f"Ill-Posed Model: Source conflict at '{comp.name}'.")

                        port.bond.causality = target_causality
                        neighbor = port.bond.target.component if port.bond.source.component == comp else port.bond.source.component
                        step1_neighbors.append(neighbor)

        self._propagate_worklist(step1_neighbors)

        # =====================================================================
        # STEP 2: Fixed Causality Type 1b (Non-Invertible / Blocks / Switched)
        # =====================================================================
        step2_neighbors: list[Component] = []
        for comp in self.components.values():
            # Targets explicit non-invertibles, signal blocks, or locked switches
            if getattr(comp, "non_invertible", False):
                for port in comp.ports.values():
                    if port.bond and port.bond.causality == Causality.UNASSIGNED:
                        fixed_causality = getattr(port, "fixed_causality", None)
                        if fixed_causality == Causality.EFFORT_AT_SOURCE:
                            port.bond.causality = (
                                Causality.EFFORT_AT_TARGET if port.bond.source == port else Causality.EFFORT_AT_SOURCE
                            )
                        elif fixed_causality == Causality.EFFORT_AT_TARGET:
                            port.bond.causality = (
                                Causality.EFFORT_AT_SOURCE if port.bond.source == port else Causality.EFFORT_AT_TARGET
                            )
                        neighbor = port.bond.target.component if port.bond.source.component == comp else port.bond.source.component
                        step2_neighbors.append(neighbor)

        self._propagate_worklist(step2_neighbors)

        # =====================================================================
        # STEP 3: Preferred Causality (Integral Causality for Storage Elements)
        # =====================================================================
        storage_types = (ComponentType.C, ComponentType.MC, ComponentType.I, ComponentType.MI, ComponentType.IC, ComponentType.MIC)
        for comp in self.components.values():
            if comp.component_type in storage_types:
                for port in comp.ports.values():
                    if port.bond:
                        # Determine if this specific port acts as C or I
                        # For mixed IC/MIC, check port-level definitions; fallback to component level
                        pref_causality = None
                        if port.port_type == PortType.C_TYPE_PORT:
                            pref_causality = Causality.EFFORT_AT_TARGET if port.bond.source == port else Causality.EFFORT_AT_SOURCE
                        elif port.port_type == PortType.I_TYPE_PORT:
                            pref_causality = Causality.EFFORT_AT_SOURCE if port.bond.source == port else Causality.EFFORT_AT_TARGET

                        if pref_causality:
                            if port.bond.causality == Causality.UNASSIGNED:
                                port.bond.causality = pref_causality
                                neighbor = port.bond.target.component if port.bond.source.component == comp else port.bond.source.component
                                self._propagate_worklist([neighbor])
                            elif port.bond.causality != pref_causality:
                                if comp not in self.derivative_causality_components:
                                    self.derivative_causality_components.append(comp)

        # =====================================================================
        # STEP 4: Arbitrary / Free Causality & Algebraic Loop Inventory
        # =====================================================================
        for comp in self.components.values():
            if comp.component_type in (ComponentType.R, ComponentType.MR):
                for port in comp.ports.values():
                    if port.bond and port.bond.causality == Causality.UNASSIGNED:
                        # Assign arbitrary effort out
                        port.bond.causality = (
                            Causality.EFFORT_AT_SOURCE if port.bond.source == port else Causality.EFFORT_AT_TARGET
                        )
                        
                        # Trace if this choice formed an algebraic loop back to itself
                        loop_bonds = self._trace_algebraic_loop(comp, port.bond)
                        if loop_bonds:
                            self.algebraic_loops.append(loop_bonds)

                        neighbor = port.bond.target.component if port.bond.source.component == comp else port.bond.source.component
                        self._propagate_worklist([neighbor])

        unassigned = [b for b in self.bonds if b.causality == Causality.UNASSIGNED]
        if unassigned:
            raise RuntimeError(f"SCAP failed: {len(unassigned)} bond(s) remained unassigned. Check ill-posed structures or disconnected loops.")

        self.system_type = self._classify_system()
        return self.system_type

    def _propagate_worklist(self, initial_components: list[Component]) -> None:
        worklist = deque(initial_components)
        while worklist:
            comp = worklist.popleft()
            self._propagate_component_causality(comp, worklist)

    def _propagate_component_causality(self, comp: Component, worklist: deque[Component]) -> bool:
        bonds = list(self._get_bonds_for_component(comp))
        assigned_bonds = [b for b in bonds if b.causality != Causality.UNASSIGNED]
        unassigned_bonds = [b for b in bonds if b.causality == Causality.UNASSIGNED]

        if not unassigned_bonds:
            return False

        changed_bonds: list[Bond] = []

        # 0-Junction & Switched 0-Junction (1 Effort IN constraint)
        if comp.component_type in (ComponentType.ZERO, ComponentType.XZERO):
            effort_in = sum(1 for b in assigned_bonds if self._get_bond_effort_direction(b, comp) == "IN")
            
            if effort_in > 1:
                raise ValueError(f"Causality Conflict: 0-Junction '{comp.name}' has {effort_in} effort inputs (max 1).")

            if effort_in == 1:
                for b in unassigned_bonds:
                    b.causality = Causality.EFFORT_AT_TARGET if b.source.component == comp else Causality.EFFORT_AT_SOURCE
                    changed_bonds.append(b)
            elif len(unassigned_bonds) == 1 and effort_in == 0:
                b = unassigned_bonds[0]
                b.causality = Causality.EFFORT_AT_SOURCE if b.source.component == comp else Causality.EFFORT_AT_TARGET
                changed_bonds.append(b)

        # 1-Junction & Switched 1-Junction (1 Effort OUT constraint)
        elif comp.component_type in (ComponentType.ONE, ComponentType.XONE):
            effort_out = sum(1 for b in assigned_bonds if self._get_bond_effort_direction(b, comp) == "OUT")

            if effort_out > 1:
                raise ValueError(f"Causality Conflict: 1-Junction '{comp.name}' has {effort_out} effort outputs / flow inputs (max 1).")

            if effort_out == 1:
                for b in unassigned_bonds:
                    b.causality = Causality.EFFORT_AT_SOURCE if b.source.component == comp else Causality.EFFORT_AT_TARGET
                    changed_bonds.append(b)
            elif len(unassigned_bonds) == 1 and effort_out == 0:
                b = unassigned_bonds[0]
                b.causality = Causality.EFFORT_AT_TARGET if b.source.component == comp else Causality.EFFORT_AT_SOURCE
                changed_bonds.append(b)

        # Transformers (TF, MTF)
        elif comp.component_type in (ComponentType.TF, ComponentType.MTF) and len(unassigned_bonds) == 1 and len(assigned_bonds) == 1:
            assigned_dir = self._get_bond_effort_direction(assigned_bonds[0], comp)
            b = unassigned_bonds[0]
            b.causality = (
                (Causality.EFFORT_AT_TARGET if b.source.component == comp else Causality.EFFORT_AT_SOURCE)
                if assigned_dir == "IN" else
                (Causality.EFFORT_AT_SOURCE if b.source.component == comp else Causality.EFFORT_AT_TARGET)
            )
            changed_bonds.append(b)

        # Gyrators (GY, MGY)
        elif comp.component_type in (ComponentType.GY, ComponentType.MGY) and len(unassigned_bonds) == 1 and len(assigned_bonds) == 1:
            assigned_dir = self._get_bond_effort_direction(assigned_bonds[0], comp)
            b = unassigned_bonds[0]
            b.causality = (
                (Causality.EFFORT_AT_SOURCE if b.source.component == comp else Causality.EFFORT_AT_TARGET)
                if assigned_dir == "IN" else
                (Causality.EFFORT_AT_TARGET if b.source.component == comp else Causality.EFFORT_AT_SOURCE)
            )
            changed_bonds.append(b)

        for b in changed_bonds:
            other_comp = b.target.component if b.source.component == comp else b.source.component
            worklist.append(other_comp)

        return len(changed_bonds) > 0

    def _trace_algebraic_loop(self, start_comp: Component, initial_bond: Bond) -> list[Bond]:
        visited_bonds = [initial_bond]
        curr = initial_bond.target.component if initial_bond.source.component == start_comp else initial_bond.source.component
        
        while curr and curr != start_comp:
            bonds = [b for b in self._get_bonds_for_component(curr) if b not in visited_bonds]
            if not bonds:
                break
            # Simply traces the first unvisited branch; robust implementations may require DFS for branching networks
            next_bond = bonds[0]
            visited_bonds.append(next_bond)
            curr = next_bond.target.component if next_bond.source.component == curr else next_bond.source.component

        return visited_bonds if curr == start_comp else []

    def _classify_system(self) -> SystemType:
        has_derivative = len(self.derivative_causality_components) > 0
        has_algebraic = len(self.algebraic_loops) > 0

        if has_derivative and has_algebraic:
            return SystemType.DAE_MIXED
        elif has_derivative:
            return SystemType.DAE_DERIVATIVE
        elif has_algebraic:
            return SystemType.DAE_ALGEBRAIC
        return SystemType.ODE

    def draw(self, filename: str = "bond_graph", format: str = "png", view: bool = True) -> graphviz.Digraph:
        """
        Renders the Bond Graph from source --> target with formal causal strokes.
        - Power flow arrow points from source to target.
        - Causal stroke is drawn at the effort-receiving end.
        """
        dot = graphviz.Digraph(name=self.name, comment="Bond Graph Visualization")
        dot.attr(rankdir="LR", nodesep="0.6", ranksep="0.8")

        # Clean textbook node styling (no bounding boxes)
        dot.attr("node", shape="plaintext", fontname="Helvetica-Bold", fontsize="14")

        # Render Component Nodes
        for comp_name, comp in self.components.items():
            if comp.component_type in (ComponentType.ONE, ComponentType.ZERO):
                label = comp.component_type.name
            else:
                label = f"{comp.component_type.name}: {comp_name}"

            dot.node(comp_name, label=label)

        # Render Bonds (Source --> Target)
        for i, bond in enumerate(self.bonds):
            src_comp = bond.source.component.name
            tgt_comp = bond.target.component.name

            # Always direct edges from source to target
            dir_style = "forward"
            power_arrow = "halfopen"

            if bond.causality == Causality.EFFORT_AT_TARGET:
                # Power arrow AND Causal stroke at target end
                arrowhead = f"tee{power_arrow}"
                arrowtail = "none"
            elif bond.causality == Causality.EFFORT_AT_SOURCE:
                # Power arrow at target end, Causal stroke at source end
                arrowhead = power_arrow
                arrowtail = "tee"
                dir_style = "both"
            else:
                # Unassigned causality (only power flow arrow)
                arrowhead = power_arrow
                arrowtail = "none"

            dot.edge(
                src_comp,
                tgt_comp,
                label=f" e{i+1}, f{i+1}",
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
 
def print_causality_table(bg: BondGraph) -> None:
    """Prints a terminal representation of bonds and causality."""
    print(f"\n--- Causality Summary: {bg.name} ---")
    print(f"{'Bond':<8} | {'Source':<12} | {'Target':<15} | {'Causality'}")
    print("-" * 60)

    for i, b in enumerate(bg.bonds, 1):
        src = b.source.component.name
        tgt = b.target.component.name

        if b.causality == Causality.EFFORT_AT_TARGET:
            direction = f"{src} |-----> {tgt}"
        elif b.causality == Causality.EFFORT_AT_SOURCE:
            direction = f"{src} <-----| {tgt}"
        else:
            direction = f"{src} ------- {tgt} (UNASSIGNED)"

        print(f"Bond {i:<3} | {src:<12} | {tgt:<15} | {direction}")

if __name__ == "__main__": 
    # Construct the system: Mass (I), Spring (C), Damper (R), Force Source (SE)
    bg = BondGraph("Mass_Spring_Damper")
    # Add components
    se = bg.add_component("SE_Force", component_type=ComponentType.SE)
    j1 = bg.add_component("J1", component_type=ComponentType.ONE)
    mass = bg.add_component("I_Mass", component_type=ComponentType.I)
    spring = bg.add_component("C_Spring", component_type=ComponentType.C)
    damper = bg.add_component("R_Damper", component_type=ComponentType.R)
    # Connect components via bonds
    bg.add_bond(se, j1)
    bg.add_bond(j1, mass)
    bg.add_bond(j1, spring)
    bg.add_bond(j1, damper)
    # Run SCAP to assign causality across all bonds
    bg.assign_causality()
    bg.draw(filename="mass_spring_damper", format="png", view=True)

    print_causality_table(bg)