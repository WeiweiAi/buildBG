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
    """Categorizes physical, junction, transducer, and control components."""
    # Borutzky, Wolfgang. Bond graph modelling of engineering systems. Vol. 103. New York: springer, 2011.
    C = auto()  # Capacitor / Spring, a storage node only contains C-type ports, the flow variable is integrated with respect to time
    I = auto()  # Inductor / Mass,  a storage node only contains I-type ports, the effort variable is integrated with respect to time
    IC = auto() # storage node contains both C-type and I-type ports
    MC = auto() # A modulated C-type storage
    MI = auto() # A modulated I-type storage
    MIC = auto() # A modulated storage node contains both C-type and I-type ports
    R = auto()  # Resistor / Damper,https://bg-rdf.org/ontologies/bondgraph-framework#Dissipator
    Re= auto() # A biochemical reaction, https://bg-rdf.org/ontologies/bondgraph-framework#Reaction
    Re_GHK= auto() # A voltage modulated biochemical reaction-Goldman-Hodgkin-Katz (GHK) ion channel
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
    # Mathematical Relations
    BLOCK = auto()  # All nodes that only have signal ports and represent mathematical relations between these signals.
    # User-defined / Custom
    CUSTOM = auto() # User-defined component   

JUNCTIONS = {ComponentType.ZERO, ComponentType.ONE, ComponentType.XZERO,  ComponentType.XONE }
class PowerVariable(Enum):
    """Identifies the effort, flow, state, and signal variables of a domain."""
    EFFORT = auto() # e.g., voltage, force, pressure
    FLOW = auto()   # e.g., current, velocity, volumetric flow rate
    QUANTITY = auto() # e.g., charge, displacement, volume
    MOMENTUM = auto() # e.g., momentum, angular momentum
    POWER = auto() # e.g., power, energy rate
    ENERGY = auto() # e.g., energy, work
    SIGNAL = auto() # A signal represents one arbitrary variable of time that may also be an effort or a flow, but not necessarily
class ConnectionType(Enum):
    """Distinguishes energy-carrying bonds from one-way signal bonds."""
    POWER_BOND = auto() # carries both effort and flow (bidirectional)
    SIGNAL_BOND = auto() # carries only signal (unidirectional)
class PortType(Enum):
    """Describes the physical role and causality behavior of a port."""
    POWER_PORT = auto() # Port for power exchange (effort and flow)
    SIGNAL_PORT = auto() # Port for signal interface (control signals)
    C_TYPE_PORT = auto() # Port for C-type storage (integrates flow to quantity), is a power port
    I_TYPE_PORT = auto() # Port for I-type storage (integrates effort to momentum), is a power port
class ConstitutiveRelationship(Enum):
    """Lists supported implicit constitutive-equation forms."""
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
    """Enumerates built-in physical domains and the abstract fallback domain."""
    ABSTRACT = auto()  # Default state: uses e, f, p, q
    ELECTRICAL = auto()
    MECHANICAL_TRANSLATIONAL = auto()
    MECHANICAL_ROTATIONAL = auto()
    HYDRAULIC = auto()
    CHEMICAL = auto()
    THERMAL = auto()
    CUSTOM = auto()
class SystemType(Enum):
    """Classifies a graph's resulting ordinary or differential-algebraic system."""
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

class DomainRegistry:
    """Central registry for domain-specific physical quantities."""

    _registry: dict[Domain | str, dict[PowerVariable | str, PhysicalQuantity]] = {
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
        }
    }

    @classmethod
    def register(cls, domain: Domain | str, variables: dict[PowerVariable | str, PhysicalQuantity]) -> None:
        """Registers a new domain."""
        if domain in cls._registry:
            raise ValueError(
                f"Domain '{domain}' is already registered."
            )  
        else:
            cls._registry[domain] = variables
    @classmethod
    def replace(cls, domain: Domain | str, variables: dict[PowerVariable | str, PhysicalQuantity]) -> None:
        """Overwrites an existing domain."""
        cls._registry[domain] = variables

    @classmethod
    def get_variables(cls, domain: 'Domain | str') -> dict[PowerVariable | str, PhysicalQuantity] | None:
        """Returns registered variable metadata for a domain, if present."""
        return cls._registry.get(domain)
@dataclass
class StateVariable:
    """Represents a time-integrated energy state of a component (q or p)."""
    variable_type: PowerVariable | str
    component: Component = field(repr=False) # Prevents Infinite Recursion Crashing
    
    @property
    def symbol(self) -> str:
        """Returns the domain-specific state symbol qualified by component name."""
        # Query the new registry
        domain_dict = DomainRegistry.get_variables(self.component.domain)
        if domain_dict and self.variable_type in domain_dict:
            base_symbol = domain_dict[self.variable_type].symbol
            return f"{base_symbol}_{self.component.name}"
        else:
            raise ValueError(f"Domain '{self.component.domain}' does not have a registered symbol for variable type '{self.variable_type}'.")
    @property
    def derivative_symbol(self) -> str:
        """Returns the time derivative of the state variable (x_dot)."""
        return f"d({self.symbol})/dt"
@dataclass
class ConstitutiveEquation:
    """Represents a single implicit relation: Phi(e, f, x, x_dot) = 0"""
    expression: Any  # Could be a string for now, or a sympy.Expr in a real solver
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
@dataclass(eq=False)
class Port:
    """Represents one typed connection point on a component."""
    label: str # the ports of each component are uniquely labelled
    component: Component = field(repr=False) # Prevents Infinite Recursion Crashing
    port_type: PortType = PortType.POWER_PORT
    domain: Domain | str = Domain.ABSTRACT
    fixed_causality: bool | None = field(default=None, repr=False)# If set, this port's causality will not be changed during causality assignment.
    causality: bool | None = field(default=None, repr=False) # True if causal stroke is at this port (receiving effort), None if unassigned;
    bond: Bond | None = field(default=None, repr=False, init=False)     

    @property
    def name(self) -> str:
        """Returns the fully qualified `<component>.<port>` identifier."""
        return f"{self.component.name}.{self.label}"    
    
    @property
    def effective_domain(self) -> Domain | str:
        """Returns the port override domain or its component's domain."""
        return self.domain if self.domain != Domain.ABSTRACT else self.component.domain
    
    def _get_symbol(self, variable_type: PowerVariable) -> str:
        """Builds a qualified variable symbol using the effective domain registry."""
        # Query the new registry
        domain_dict = DomainRegistry.get_variables(self.effective_domain)
        if domain_dict and variable_type in domain_dict:
            return f"{domain_dict[variable_type].symbol}_{self.name}"
        raise ValueError(f"Domain '{self.component.domain}' does not have a registered symbol for variable type '{variable_type}'.")

    def _attach_bond(self, bond: Bond) -> None:
        """Attaches a bond to this port."""
        self.bond = bond
        self.component.hold_port(self) # Mark the port as held when a bond is attached
        self.component.bonds.add(bond) # Add the bond to the component's bond set

    def _detach_bond(self, bond: Bond) -> None:
        """Detaches the bond from this port."""
        if self.bond is not bond:
            raise ValueError( f"Port '{self.name}' is not connected to the specified bond." )
        self.component.bonds.discard(self.bond) # Remove the bond from the component's bond set
        self.bond = None
        self.component.release_port(self) # Mark the port as free when a bond is detached        

    @property
    def effort(self) -> str:
        """Returns this port's effort-variable symbol."""
        return self._get_symbol(PowerVariable.EFFORT)
    @property
    def flow(self) -> str:
        """Returns this port's flow-variable symbol."""
        return self._get_symbol(PowerVariable.FLOW)
    @property
    def quantity(self) -> str:
        """Returns this port's quantity-variable symbol."""
        return self._get_symbol(PowerVariable.QUANTITY)
    @property
    def momentum(self) -> str:
        """Returns this port's momentum-variable symbol."""
        return self._get_symbol(PowerVariable.MOMENTUM)

    @property
    def signal(self) -> str:
        """Returns this port's signal-variable symbol."""
        return self._get_symbol(PowerVariable.SIGNAL)

@dataclass(eq=False)
class Bond:
    """Connects two ports and owns their shared causality assignment."""
    source: Port
    target: Port
    connection_type: ConnectionType = ConnectionType.POWER_BOND

    def __post_init__(self) -> None:
        """Validates endpoints and atomically attaches the bond to both ports."""

        self.source._attach_bond(self)
        self.target._attach_bond(self)

    def get_other_port(self, port: Port) -> Port:
        """Returns the opposite port of the bond given one endpoint."""
        if port is self.source:
            return self.target
        elif port is self.target:
            return self.source
        else:
            raise ValueError(
                f"Port '{port.name}' is not connected to this bond."
            )
    def get_other_component(self, component: Component) -> Component:
        """Returns the opposite component of the bond given one endpoint."""
        if component is self.source.component:
            return self.target.component
        elif component is self.target.component:
            return self.source.component
        else:
            raise ValueError(
                f"Component '{component.name}' is not connected to this bond."
            )
    def validate_causality(self) -> None:
        source = self.source.causality
        target = self.target.causality
    
        if source is None and target is None:
            return
    
        if source is None or target is None:
            raise ValueError(
                f"Bond '{self.name}' has partially assigned causality."
            )
    
        if source == target:
            raise ValueError(
                f"Bond '{self.name}' has conflicting causality: "
                "both endpoints have the same causality."
            )
    def has_causality_conflict(self) -> bool | None:
        if self.source.causality is None or self.target.causality is None:
            return None # Causality is unassigned for at least one port
        return self.source.causality == self.target.causality

    def assign_causality(self,port: Port, target_causality: bool) -> bool:
        """Assigns causality to the bond's source and target ports, given a target causality for the specified port."""
        if check := self.has_causality_conflict():
            raise ValueError(
                f"Cannot assign causality: source '{self.source.name}' and target '{self.target.name}' have conflicting causality assignments or unassigned ports."
            )
        other_port = self.get_other_port(port)

        if port.causality is not None:
            if port.causality != target_causality:
                raise ValueError(f"Cannot assign causality: port '{port.name}' already has a conflicting causality assignment.")

        if port.fixed_causality is not None:
            if port.fixed_causality != target_causality:
                raise ValueError(f"Cannot assign causality: port '{port.name}' has a fixed causality that conflicts with the target causality.")
        other_causality = not target_causality

        if other_port.causality is not None:
            if other_port.causality != other_causality:
                raise ValueError(f"Cannot assign causality: port '{other_port.name}' already has a conflicting causality assignment.")

        if other_port.fixed_causality is not None:
            if other_port.fixed_causality != other_causality:
                raise ValueError(f"Cannot assign causality: port '{other_port.name}' has a fixed causality that conflicts with the target causality.")

        # Commit only after all validation succeeds.
        port.causality = target_causality
        other_port.causality = other_causality

        return True # Successfully assigned causality to the specified port; the other port's causality will be the opposite.

    @property
    def name(self) -> str:
        """Returns the source-to-target bond identifier."""
        return f"{self.source.name}--{self.target.name}"

    def disconnect(self) -> None:
        """Safely severs the bidirectional link between the bond and its ports."""
        self.source._detach_bond(self)
        self.target._detach_bond(self)

    def validate_connection(self) -> None:
        """Validates the bond's source and target ports."""
        Bond.validate(self.source, self.target, self.connection_type)
    
    @staticmethod
    def validate(
        source: Port,
        target: Port,
        connection_type: ConnectionType = ConnectionType.POWER_BOND,
    ) -> None:

        if source is target:
            raise ValueError(
                f"Cannot create a bond from port '{source.name}' to itself."
            )

        if source.component is target.component:
            raise ValueError(
                "Cannot create a bond between ports on the same component."
            )

        if source.bond is not None:
            raise ValueError(
                f"Source port {source.name} is already connected to a bond."
            )

        if target.bond is not None:
            raise ValueError(
                f"Target port {target.name} is already connected to a bond."
            )

        if connection_type == ConnectionType.POWER_BOND:
            valid_power_types = (
                PortType.POWER_PORT,
                PortType.C_TYPE_PORT,
                PortType.I_TYPE_PORT,
            )

            if source.port_type not in valid_power_types:
                raise ValueError(
                    f"Source port {source.name} is not a valid power port."
                )

            if target.port_type not in valid_power_types:
                raise ValueError(
                    f"Target port {target.name} is not a valid power port."
                )

            if source.effective_domain != target.effective_domain:
                raise ValueError(
                    f"Domain mismatch: {source.name} "
                    f"({source.effective_domain}) cannot be connected to "
                    f"{target.name} ({target.effective_domain})."
                )

        elif connection_type == ConnectionType.SIGNAL_BOND:
            if source.port_type is not PortType.SIGNAL_PORT:
                raise ValueError(
                    f"Source port {source.name} is not a valid signal port."
                )

            if target.port_type is not PortType.SIGNAL_PORT:
                raise ValueError(
                    f"Target port {target.name} is not a valid signal port."
                )

        else:
            raise ValueError(
                f"Unsupported connection type: {connection_type!r}."
            )

    @property
    def effort(self) -> str:
        """Returns the effort symbol supplied by the causality assignment."""
        self.validate_connection()  # Ensure the bond is valid before accessing effort
        if self.source.causality : # Causal stroke at the source port
            return self.target.effort
        elif self.target.causality: # Causal stroke at the target port
            return self.source.effort
        else:
            if self.connection_type == ConnectionType.SIGNAL_BOND:
                # For signal bonds, we can default to source effort if causality is unassigned
                return self.source.effort
            else:
                raise ValueError("Causality is not assigned.")
        
    @property
    def flow(self) -> str:
        """Returns the flow symbol supplied by the causality assignment."""
        self.validate_connection()  # Ensure the bond is valid before accessing flow
        if self.source.causality:
            return self.source.flow
        elif self.target.causality:
            return self.target.flow
        else:
            if self.connection_type == ConnectionType.SIGNAL_BOND:
                # For signal bonds, we can default to source flow if causality is unassigned
                return self.source.flow
            else:
                raise ValueError("Causality is not assigned.")
    @property
    def signal(self) -> str:
        """Returns the source signal for a signal bond."""
        if self.connection_type != ConnectionType.SIGNAL_BOND:
            raise ValueError("Signal property is only valid for signal bonds.")
        return self.source.signal

@dataclass(eq=False)
class Component:
    """Models a bond-graph element, its ports, parameters, states, and equations."""
    name: str
    component_type: ComponentType | str = ComponentType.CUSTOM
    domain: Domain | str = Domain.ABSTRACT
    non_invertible: bool = False # If True, the component has any constitutive relationship that cannot be algebraically inverted to solve for either effort or flow.
    # Optional parameters strictly for ComponentType.CUSTOM
    custom_power_ports: list[bool | None] = field(default_factory=list) # List of Causality enums for each custom power port
    custom_signal_ports: int = 0
    ports: dict[str, Port] = field(default_factory=dict, repr=False, init=False)
    _available_ports: deque[Port] = field(default_factory=deque, repr=False,   init=False) # track which ports are available for new bonds  
    _next_port_number: int = field(default=1, repr=False, init=False) # only used for junctions, to auto-label new ports
    bonds: set[Bond] = field(default_factory=set, repr=False, init=False) 
    # consider the following later  
    parameters: dict[str, Any]  = field(default_factory=dict, repr=False, init=False)    
    states: dict[PowerVariable, StateVariable]  = field(default_factory=dict, repr=False, init=False)
    equations: list[ConstitutiveEquation] = field(default_factory=list, repr=False, init=False)
    # Optional callback to override the default linear equations
    equation_generator: Callable[[Component], list[ConstitutiveEquation]] | None = None
        
    def __post_init__(self) -> None:
        """Creates default ports based on the component type and requested count."""
        # 1-Port Elements
        if self.component_type in (ComponentType.R, ComponentType.SE, ComponentType.SF ):
            p1=self._add_port("p1")
            self._available_ports.append(p1) # For 1-port elements, the single port is always available for bonding
        elif self.component_type == ComponentType.C:
            p1=self._add_port("p1", port_type=PortType.C_TYPE_PORT)
            self._available_ports.append(p1)
        elif self.component_type == ComponentType.I:
            p1=self._add_port("p1", port_type=PortType.I_TYPE_PORT)
            self._available_ports.append(p1)
        # 2-Port Elements
        elif self.component_type in (ComponentType.TF, ComponentType.GY):
            p1=self._add_port("p1")
            p2=self._add_port("p2")
            self._available_ports.extend([p1, p2])
        elif self.component_type in (ComponentType.MSE, ComponentType.MSF):
            p1=self._add_port("p1")
            p2=self._add_port("mod", port_type=PortType.SIGNAL_PORT)
            self._available_ports.extend([p1, p2])
        elif self.component_type == ComponentType.IC:
            p1=self._add_port("p1", port_type=PortType.I_TYPE_PORT)
            p2=self._add_port("p2", port_type=PortType.C_TYPE_PORT)
            self._available_ports.extend([p1, p2])
        elif self.component_type == ComponentType.MC:
            p1=self._add_port("p1", port_type=PortType.C_TYPE_PORT)
            p2=self._add_port("mod", port_type=PortType.SIGNAL_PORT)
            self._available_ports.extend([p1, p2])
        elif self.component_type == ComponentType.MI:
            p1=self._add_port("p1", port_type=PortType.I_TYPE_PORT)
            p2=self._add_port("mod", port_type=PortType.SIGNAL_PORT)
            self._available_ports.extend([p1, p2])
        # Reaction Elements
        elif self.component_type == ComponentType.Re:
            p1=self._add_port("p1",fixed_causality=True)
            p2=self._add_port("p2", fixed_causality=True)
            self._available_ports.extend([p1, p2])
            self.non_invertible = True # Reactions are generally non-invertible due to their nonlinear constitutive relationships
        # 3-Port Elements
        elif self.component_type == ComponentType.MIC:
            p1=self._add_port("p1", port_type=PortType.I_TYPE_PORT)
            p2=self._add_port("p2", port_type=PortType.C_TYPE_PORT)
            p3=self._add_port("mod", port_type=PortType.SIGNAL_PORT)
            self._available_ports.extend([p1, p2, p3])
        elif self.component_type == ComponentType.Re_GHK:
            p1=self._add_port("p1", fixed_causality=True)
            p2=self._add_port("p2", fixed_causality=True)
            p3=self._add_port("mod", port_type=PortType.SIGNAL_PORT)
            self._available_ports.extend([p1, p2, p3])
            self.non_invertible = True # Modulated storage elements are generally non-invertible due to their nonlinear constitutive relationships
        elif self.component_type in (ComponentType.MTF, ComponentType.MGY):
            p1=self._add_port("p1")
            p2=self._add_port("p2")
            p3=self._add_port("mod", port_type=PortType.SIGNAL_PORT)
            self._available_ports.extend([p1, p2, p3])
        elif self.component_type in (ComponentType.ZERO, ComponentType.ONE, ComponentType.XZERO, ComponentType.XONE):
            pass # Junctions dynamically allocate ports as needed; no default ports are created.
        elif self.component_type == ComponentType.BLOCK:
            for i in range(1, self.custom_signal_ports + 1):
                p = self._add_port(f"s{i}", port_type=PortType.SIGNAL_PORT)
                self._available_ports.append(p)
        else: # component_type == ComponentType.CUSTOM or any other unrecognized type 
            # For custom components, create the specified number of power and signal ports
            # check that custom_power_ports is a list of Causality or None
            if not isinstance(self.custom_power_ports, list) or not all(isinstance(c, (bool, type(None))) for c in self.custom_power_ports):
                raise ValueError("custom_power_ports must be a list of bool or None.")
            # if there is any fixed causality in custom_power_ports, then the component is non-invertible
            if any(c is not None for c in self.custom_power_ports):
                self.non_invertible = True
            for i in range(1, len(self.custom_power_ports) + 1):
                p = self._add_port(f"p{i}", port_type=PortType.POWER_PORT, fixed_causality=self.custom_power_ports[i - 1])
                self._available_ports.append(p)
            for i in range(1, self.custom_signal_ports + 1):
                p = self._add_port(f"s{i}", port_type=PortType.SIGNAL_PORT)
                self._available_ports.append(p)

    @property
    def port_count(self) -> int:
        """Dynamically always returns the true number of ports."""
        return len(self.ports)

    def _add_port(self, label: str, **kwargs) -> Port:
        if label in self.ports: # uniqueness check for port labels
            raise ValueError(f"Port '{label}' already exists on component '{self.name}'.")       
        new_port = Port(label=label, component=self, **kwargs)
        self.ports[label] = new_port
        return new_port

    def release_port(self, port: Port) -> None:
        if port.bond is not None:
            raise ValueError(
                f"Cannot release connected port '{port.name}'."
            )
        if port.component is not self:
            raise ValueError(
                f"Port '{port.name}' does not belong to component '{self.name}'."
            )    
        if port not in self._available_ports:
            self._available_ports.append(port)
        else:
            warnings.warn(f"Port '{port.name}' is already marked as available."); 

    def hold_port(self, port: Port) -> None:
        if port.bond is None:
            raise ValueError(
                f"Cannot hold unconnected port '{port.name}'."
            )
        if port.component is not self:
            raise ValueError(
                f"Port '{port.name}' does not belong to component '{self.name}'."
            )    
        if port in self._available_ports:
            self._available_ports.remove(port)
        else:
            warnings.warn(f"Port '{port.name}' is already marked as held.");
    
    def get_or_create_port(self) -> Port:
        # This method is only relevant for junctions (0, 1, X0, X1). It creates a new one port.
        if self.component_type in (
            ComponentType.ZERO,
            ComponentType.ONE,
            ComponentType.XZERO,
            ComponentType.XONE,
        ):
            if self._available_ports:
                return self._available_ports[0]  # Return the first available free port
            else: # Create a new one.
                label = f"p{self._next_port_number}" 
                self._next_port_number += 1
                port = self._add_port(label)
                self._available_ports.append(port)  # Mark the new port as available for bonding              
                return port
        else:
            raise ValueError(
                f"Component '{self.name}' of type '{self.component_type}' does not support dynamic port allocation."
            )   

    def get_available_ports(self) -> list[Port]:
        """Returns a list of currently unconnected ports."""
        # This method is only relevant for non-junction components. For junctions, use `allocate_port()` to get a free port or create a new one.
        return list(self._available_ports)
class BondGraph:
    """Owns a connected set of components and assigns bond causalities."""

    def __init__(self, name: str = "bond_graph") -> None:
        """Initializes an empty named graph and causality diagnostics."""
        self.name = name
        self.components: dict[str, Component] = {} # Mapping of component names to Component objects
        # Insertion-ordered mapping gives O(1) bond membership/deletion
        # while retaining deterministic iteration order.
        self._bonds: dict[Bond, None] = {}

        # Extended Diagnostic State
        self.derivative_causality_components: list[Component] = []
        self.algebraic_loops: list[list[Bond]] = []
        self.system_type: SystemType = SystemType.ODE

    def _resolve_string(self, arg: str) -> Port | Component |None:
        """Resolves a String input into a valid Port object """
        if "." in arg:
            comp_name, port_label = arg.split(".", 1)
            comp = self.components.get(comp_name)
            if comp and port_label in comp.ports:
                return comp.ports[port_label]
            else:
                warnings.warn(f"Port '{port_label}' not found on component '{comp_name}'.")
                return None               
        else:
            comp = self.components.get(arg)
            if comp:
                return comp
            else:
                warnings.warn(f"Component '{arg}' not found in bond graph.")
                return None  
    
    @property
    def bonds(self):
        """Insertion-ordered, set-like view of all bonds."""
        return self._bonds.keys()

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
                kind="bond"
            )
            
        return G

    def add_component(self, component: Component | str, **kwargs) -> Component:
        """Adds a component to internal tracking."""
        if isinstance(component, str):
            if component in self.components:
                raise ValueError(f"Component '{component}' already exists.")
            else:
                comp_obj = Component(name=component, **kwargs)
        else:
            if component.name in self.components:
                raise ValueError(f"Component '{component.name}' already exists.")
            else:
                comp_obj = component
        # Add the component to the graph's registry
        self.components[comp_obj.name] = comp_obj
        return comp_obj

    def _resolve_endpoint(self, arg: Port | Component | str) -> Port | None:
        """Resolves a Port, Component, or String input into a valid Port object before adding a bond."""
        if isinstance(arg, Port):
            return arg
        elif isinstance(arg, Component):
            available_ports = arg.get_available_ports()
            if len(available_ports) == 1:
                return available_ports[0]  # Return the first available free port
            else:
                # For junctions, allocate a free port or create a new one
                if arg.component_type in (ComponentType.ZERO, ComponentType.ONE, ComponentType.XZERO, ComponentType.XONE):
                    return arg.get_or_create_port()  # Dynamically allocate a new port if none are free
                else:
                    warnings.warn(f"Component '{arg.name}' has {len(available_ports)} free ports, please specify which one to use.")
                    return None
        elif isinstance(arg, str):
            resolved = self._resolve_string(arg)
            if isinstance(resolved, Port):
                return resolved
            elif isinstance(resolved, Component):
                return self._resolve_endpoint(resolved)  # Recursively resolve the component to a port
            else:
                warnings.warn(f"Could not resolve '{arg}' to a valid port or component.")
                return None
        else:
            warnings.warn(f"Invalid argument type: {type(arg)}. Expected Port, Component, or str.")
            return None
   
    def add_bond(self, source: Component | Port | str, target: Component | Port | str, **kwargs) -> Bond | None:
        """Creates a bond between two endpoints."""
        src_port = self._resolve_endpoint(source)
        tgt_port = self._resolve_endpoint(target)     
        if src_port is not None and tgt_port is not None:
            try:
                Bond.validate(src_port, tgt_port, kwargs.get('connection_type', ConnectionType.POWER_BOND))
            except ValueError as e:
                warnings.warn(f"Failed to create bond: {e}")
                return None            
            bond = Bond(source=src_port, target=tgt_port, **kwargs)
            self._bonds[bond] = None   
            return bond          
        else:
            warnings.warn("Could not resolve both source and target ports for bond creation.")
            return None
       
    def _get_bonds_for_component(self, comp_arg: Component | str) -> set[Bond]:
        """Returns all bonds connected to a given component."""
        comp_name = comp_arg.name if isinstance(comp_arg, Component) else comp_arg
        comp = self.components.get(comp_name)
        if not comp:
            warnings.warn(f"Component '{comp_name}' not found in bond graph.")
            return set()
        else:
            return comp.bonds
    
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
            if (isinstance(arg1, Port) or isinstance(arg1, str)) and \
               (isinstance(arg2, Port) or isinstance(arg2, str)):
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
        for bond in bonds_to_delete:
            if bond in self._bonds:            
               # 1. Ask the bond to unhook itself from its ports
               bond.disconnect()            
               # 2. Remove it from the central graph registry
               del self._bonds[bond]
               deleted_count += 1  
            else:
                warnings.warn(f"Bond '{bond.name}' not found in the graph. It may have already been deleted.")             

        return deleted_count

    def delete_component(self, comp_arg: Component | str) -> None:
        """Removes a component, its ports, and attached bonds from the graph."""
        comp_name = comp_arg.name if isinstance(comp_arg, Component) else comp_arg
        comp = self.components.get(comp_name)
        if not comp:
            warnings.warn(f"Component '{comp_name}' not found in bond graph.")
            return

        # delete_bond() mutates Component.bonds, so iterate over a snapshot.
        for bond in list(self._get_bonds_for_component(comp)):
            self.delete_bond(bond)

        del self.components[comp.name]

    def _get_bond_effort_direction(self, bond: Bond, component: Component) -> str | None:
        """Determines if effort is flowing 'IN' to or 'OUT' of the given component via this bond."""
        if bond.source.causality is None or bond.target.causality is None:
            return None # Causality not yet assigned for this bond
            
        if bond.source.component == component:
            return "IN" if bond.source.causality  else "OUT"
        elif bond.target.component == component:
            return "IN" if bond.target.causality  else "OUT"
        return None

    def assign_causality(self) -> SystemType:
        """Executes the Generalized Extended SCAP framework."""

        for comp in self.components.values():
            for port in comp.ports.values():
                port.causality = port.fixed_causality   

        self.derivative_causality_components = []
        self.algebraic_loops = []

        # =====================================================================
        # STEP 1: Fixed Causality Type 1a (Independent & Modulated Sources)
        # =====================================================================
        step123_neighbors: list[Component] = []
        for comp in self.components.values():
            if comp.component_type in (ComponentType.SE, ComponentType.MSE, ComponentType.SF, ComponentType.MSF):
                for port in comp.ports.values():
                    if port.bond: # active bond
                        is_effort_source = comp.component_type in (ComponentType.SE, ComponentType.MSE)
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
        for comp in self.components.values():
            # Targets explicit non-invertibles, signal blocks, or locked switches
            if getattr(comp, "non_invertible", False):
                for port in comp.ports.values():
                    if port.bond :
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
        for comp in self.components.values():
            if comp.component_type in storage_types:
                for port in comp.ports.values():
                    if port.bond:
                        # Determine if this specific port acts as C or I
                        # For mixed IC/MIC, check port-level definitions; fallback to component level
                        pref_causality = None
                        if port.port_type == PortType.C_TYPE_PORT:
                            pref_causality = False
                        elif port.port_type == PortType.I_TYPE_PORT:
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
        for comp in self.components.values():
            if comp.component_type in (ComponentType.R, ComponentType.MR):
                for port in comp.ports.values():
                    if port.bond :
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

        unassigned = [b for b in self.bonds if b.source.causality is None or b.target.causality is None]
        if unassigned:
            print(f"Unassigned bonds: {[b.name for b in unassigned]}")
            raise RuntimeError(f"SCAP failed: {len(unassigned)} bond(s) remained unassigned. Check ill-posed structures or disconnected loops.")
        for bond in self.bonds:
            bond.validate_causality() # Ensure all bonds are valid after causality assignment
        self.system_type = self._classify_system()
        return self.system_type

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

    def _propagate_component_causality(
        self,
        comp: Component,
        worklist: deque[Component],
        queued: set[Component],
    ) -> bool:
        """Applies the local junction or transducer causality rules for one component."""
        
        assigned_bonds: list[Bond] = []
        unassigned_bonds: list[Bond] = []
        for bond in self._get_bonds_for_component(comp):
            if bond.source.causality is None or bond.target.causality is None:
                unassigned_bonds.append(bond)
            else:
                assigned_bonds.append(bond)

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
                    b.assign_causality(b.source, False) if b.source.component == comp else b.assign_causality(b.source, True)
                    changed_bonds.append(b)
            elif len(unassigned_bonds) == 1 and effort_in == 0:
                b = unassigned_bonds[0]
                b.assign_causality(b.source, False) if b.source.component == comp else b.assign_causality(b.source, True)
                changed_bonds.append(b)

        # 1-Junction & Switched 1-Junction (1 Effort OUT constraint)
        elif comp.component_type in (ComponentType.ONE, ComponentType.XONE):
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
        elif comp.component_type in (ComponentType.TF, ComponentType.MTF) and len(unassigned_bonds) == 1 and len(assigned_bonds) == 1:
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
        elif comp.component_type in (ComponentType.GY, ComponentType.MGY) and len(unassigned_bonds) == 1 and len(assigned_bonds) == 1:
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
                (b for b in self._get_bonds_for_component(curr) if b not in visited_bonds),
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
        

        for comp in self.components.values():
            if comp.component_type in JUNCTIONS:
                continue

            for port in comp.ports.values():
                if port.bond is None:
                    continue

                if port.causality is True:
                    # Component receives effort and provides flow.
                    variable_expr[port]["flow"] = port.flow

                elif port.causality is False:
                    # Component receives flow and provides effort.
                    variable_expr[port]["effort"] = port.effort

                else:
                    raise RuntimeError(
                        f"Port '{port.name}' has no assigned causality."
                    )
    def _propagate_across_bond(
            self,
            port: Port,
            variable: str,
            expression: str,
            variable_expr: dict,
        ) -> bool:

            bond = port.bond

            if bond is None:
                raise RuntimeError(
                    f"Port '{port.name}' is not connected to a bond."
                )

            other_port = bond.get_other_port(port)

            existing = variable_expr[other_port][variable]

            if existing is None:
                variable_expr[other_port][variable] = expression
                return True

            if existing != expression:
                raise ValueError(
                    f"Conflicting {variable} expressions at "
                    f"port '{other_port.name}': "
                    f"'{existing}' vs '{expression}'."
                )

            return False 
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
                return True

            if current != expression:
                raise ValueError(
                    f"Conflicting {variable} expressions for "
                    f"port '{port.name}': "
                    f"'{current}' versus '{expression}'."
                )

            return False
    def _propagate_junction(
            self,
            junction: Component,
            rule: JunctionRule,
            variable_expr: dict,
        ) -> JunctionPropagation | None:

            power_bonds = [
                bond
                for bond in junction.bonds
                if bond.connection_type == ConnectionType.POWER_BOND
            ]

            if len(power_bonds) < 2:
                return None

            junction_ports = [
                self._get_port_for_component(bond, junction)
                for bond in power_bonds
            ]

            source_ports = [
                port
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

            expression = variable_expr[source_other_port][
                rule.propagated_variable
            ]

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
    def _generate_component_input_equations(
            self,
            variable_expr,
        ) -> list[ConstitutiveEquation]:

            equations = []

            junction_types = {
                ComponentType.ZERO,
                ComponentType.ONE,
                ComponentType.XZERO,
                ComponentType.XONE,
            }

            for comp in self.components.values():
            
                if comp.component_type in junction_types:
                    continue
                
                for port in comp.ports.values():
                
                    if port.bond is None:
                        continue
                    
                    if port.causality is True:
                        expr = variable_expr[port]["effort"]

                        if expr is None:
                           continue

                        equations.append(
                            ConstitutiveEquation(
                                f"{port.effort} = ({expr})",
                                f"Network effort equation for {port.name}",
                            )
                        )

                    elif port.causality is False:
                        expr = variable_expr[port]["flow"]

                        if expr is None:
                            continue

                        equations.append(
                            ConstitutiveEquation(
                                f"{port.flow} = ({expr})",
                                f"Network flow equation for {port.name}",
                            )
                        )

            return equations 
    def _junction_conservation(
            self,
            junction: Component,
            rule: JunctionRule,
            state: JunctionPropagation,
        ) -> ConstitutiveEquation | None:
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

                expression = getattr(
                    other_port,
                    rule.conserved_variable,
                )

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
                    f"{prefix}({expression})"
                )

            rhs = "".join(rhs_terms)

            if not rhs:
                rhs = "0"

            return ConstitutiveEquation(
                f"({dependent_expression}) = {rhs}",
                f"{junction.component_type.name}-junction "
                f"'{junction.name}': "
                f"{rule.conserved_variable} conservation",
            )
    def _resolve_variable_expression(
        self,
        port: Port,
        variable: str,
        variable_expr: dict,
    ) -> str:

        expression = variable_expr[port][variable]

        if expression is not None:
            return expression

        if port.component.component_type not in (
            ComponentType.ZERO,
            ComponentType.ONE,
            ComponentType.XZERO,
            ComponentType.XONE,
        ):
            # For a non-junction, its own symbol is allowed.
            return getattr(port, variable)

        raise RuntimeError(
            f"Could not resolve {variable} for junction port "
            f"'{port.name}' without using a junction variable."
        )
    def _generate_junction_equations(
            self,
            junction: Component,
            variable_expr: dict,
        ) -> list[ConstitutiveEquation]:

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
    def _propagate_junction_variables(
            self,
            variable_expr: dict,
        ) -> dict[Component, JunctionPropagation]:
            """
            Propagate known junction variables until no further propagation
            is possible.

            Returns the propagation state for each junction that was
            successfully resolved.
            """

            junctions = [
                comp
                for comp in self.components.values()
                if comp.component_type in (
                    ComponentType.ZERO,
                    ComponentType.ONE,
                    ComponentType.XZERO,
                    ComponentType.XONE,
                )
            ]

            propagation_states: dict[Component, JunctionPropagation] = {}

            changed = True

            while changed:
                changed = False

                for junction in junctions:
                
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
                        continue
                    
                    state = self._propagate_junction(
                        junction,
                        rule,
                        variable_expr,
                    )

                    if state is not None:
                        propagation_states[junction] = state

                        # The generic propagation function should return
                        # whether any new variable expression was established.
                        if state.changed:
                            changed = True

            return propagation_states
    def generate_network_equations(self) -> list[ConstitutiveEquation]:
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
            for bond in self.bonds:
                if bond.source.causality is None or bond.target.causality is None:
                    raise RuntimeError(
                        f"Cannot generate equations: bond '{bond.name}' "
                        "has unassigned causality."
                    )

                if bond.source.causality == bond.target.causality:
                    raise RuntimeError(
                        f"Cannot generate equations: bond '{bond.name}' "
                        "has invalid causality; both endpoints have the "
                        f"same value ({bond.source.causality})."
                    )

            # ------------------------------------------------------------------
            # 1. Create the propagation map.
            #
            # None means:
            #   "this variable has not been determined by network propagation."
            #
            # We deliberately do NOT use "" because an empty string is not
            # semantically different from a missing expression.
            # ------------------------------------------------------------------
            variable_expr: dict[Port, dict[str, str | None]] = {
                port: {
                    "effort": None,
                    "flow": None,
                }
                for comp in self.components.values()
                for port in comp.ports.values()
                if port.bond is not None
            }

            # ------------------------------------------------------------------
            # 2. Initialize component-side known variables.
            #
            # True  -> receives effort -> flow is component-side known
            # False -> provides effort -> effort is component-side known
            # ------------------------------------------------------------------
           
            for comp in self.components.values():
                if comp.component_type in JUNCTIONS:
                    continue
                
                for port in comp.ports.values():
                    if port.bond is None:
                        continue
                    
                    if port.causality is True:
                        # Effort will be determined by the network.
                        # Flow is supplied by the component.
                        variable_expr[port]["flow"] = port.flow

                    elif port.causality is False:
                        # Flow will be determined by the network.
                        # Effort is supplied by the component.
                        variable_expr[port]["effort"] = port.effort

                    else:
                        raise RuntimeError(
                            f"Port '{port.name}' has unassigned causality."
                        )

            # ------------------------------------------------------------------
            # 3. Propagate junction variables.
            #
            # propagation_states[junction] contains the dependent junction
            # port and its corresponding non-junction port, so the conservation
            # equation can reuse the result rather than searching again.
            # ------------------------------------------------------------------
            propagation_states = self._propagate_junction_variables(
                variable_expr
            )

            # ------------------------------------------------------------------
            # 4. Generate network input equations for non-junction components.
            #
            # True  -> effort is network supplied
            # False -> flow is network supplied
            #
            # If propagation did not determine the required variable, do NOT
            # generate an invalid equation such as:
            #
            #     e_I_Mass.p1 - ()
            #
            # That variable will instead appear in a conservation equation.
            # ------------------------------------------------------------------
            equations: list[ConstitutiveEquation] = []

            for comp in self.components.values():
                if comp.component_type in JUNCTIONS:
                    continue
                
                for port in comp.ports.values():
                    if port.bond is None:
                        continue
                    
                    if port.causality is True:
                        # Port receives effort.
                        # Therefore effort is supplied by the network.
                        expression = variable_expr[port]["effort"]

                        if expression is not None:
                            equations.append(
                                ConstitutiveEquation(
                                    f"({port.effort}) = ({expression})",
                                    f"Network effort equation for {port.name}",
                                )
                            )

                    elif port.causality is False:
                        # Port receives flow.
                        # Therefore flow is supplied by the network.
                        expression = variable_expr[port]["flow"]

                        if expression is not None:
                            equations.append(
                                ConstitutiveEquation(
                                    f"({port.flow}) = ({expression})",
                                    f"Network flow equation for {port.name}",
                                )
                            )

            # ------------------------------------------------------------------
            # 5. Generate junction conservation equations.
            # ------------------------------------------------------------------
            for junction, propagation in propagation_states.items():
            
                if junction.component_type in (
                    ComponentType.ZERO,
                    ComponentType.XZERO,
                ):
                    equation = self._junction_conservation(
                        junction=junction,
                        rule=ZERO_RULE,
                        state=propagation,
                    )

                elif junction.component_type in (
                    ComponentType.ONE,
                    ComponentType.XONE,
                ):
                    equation = self._junction_conservation(
                        junction=junction,
                        rule=ONE_RULE,
                        state=propagation,
                    )

                else:
                    continue
                
                if equation is not None:
                    equations.append(equation)

            return equations                         
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
            label=''
            if comp.component_type in (ComponentType.ONE, ComponentType.ZERO, ComponentType.XONE, ComponentType.XZERO):
                if comp.component_type == ComponentType.ONE:
                    label = "1"
                elif comp.component_type == ComponentType.ZERO:
                    label = "0"
                elif comp.component_type == ComponentType.XONE:
                    label = "X1"
                elif comp.component_type == ComponentType.XZERO:
                    label = "X0"
            else:
                if isinstance(comp.component_type, ComponentType):
                    label = f"{comp.component_type.name}: {comp_name}"
                else:
                    label = f"{comp.component_type}: {comp_name}"

            dot.node(comp_name, label=label)

        # Render Bonds (Source --> Target)
        for i, bond in enumerate(self.bonds):
            src_comp = bond.source.component.name
            tgt_comp = bond.target.component.name

            # Always direct edges from source to target
            dir_style = "forward"
            power_arrow = "halfopen"

            if bond.target.causality == True:
                # Power arrow AND Causal stroke at target end
                arrowhead = f"tee{power_arrow}"
                arrowtail = "none"
            elif bond.source.causality == True:
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

        if b.target.causality == True:
            direction = f"{src} |-----> {tgt}"
        elif b.source.causality == True:
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
    bg_equations = bg.generate_network_equations()
    for eq in bg_equations:
        print(eq)