from __future__ import annotations
from enum import Enum, auto
from dataclasses import asdict, dataclass, field, is_dataclass
from collections import deque
import warnings
import json

# Define basic data structures for bond graph modeling, including components, ports, bonds, and the overall bond graph.
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
TRANSDUCERS = {ComponentType.TF, ComponentType.GY, ComponentType.MTF, ComponentType.MGY }

class ConnectionType(Enum):
    """Distinguishes energy-carrying bonds from one-way signal bonds."""
    POWER_BOND = auto() # carries both effort and flow (bidirectional)
    SIGNAL_BOND = auto() # carries only signal (unidirectional)
class PortType(Enum):
    """Describes the physical role and causality behavior of a port."""
    POWER_PORT = auto() # Port for power exchange (effort and flow)
    SIGNAL_PORT = auto() # Port for signal interface (control signals)
class StorageType(Enum):
    C_TYPE = auto() # Port for C-type storage (integrates flow to quantity), is a power port
    I_TYPE = auto() # Port for I-type storage (integrates effort to momentum), is a power port
class Domain(Enum):
    """Enumerates built-in physical domains and the abstract fallback domain."""
    ABSTRACT = auto()  # Default state: uses e, f, p, q
    ELECTRICAL = auto()
    MECHANICAL_TRANSLATIONAL = auto()
    MECHANICAL_ROTATIONAL = auto()
    HYDRAULIC = auto()
    CHEMICAL = auto()
    ELECTROCHEMICAL = auto()
    THERMAL = auto()
    CUSTOM = auto()
@dataclass(frozen=True)
class PhysicalQuantity:
    """Metadata for a domain-specific power variable."""
    description: str
    symbol: str
    units: str
    value: float | None = None  # Optional numerical value for the physical quantity

class BGVariable(Enum):
    """Identifies the effort, flow, state, and signal variables of a domain."""
    EFFORT = auto() # e.g., voltage, force, pressure
    FLOW = auto()   # e.g., current, velocity, volumetric flow rate
    QUANTITY = auto() # e.g., charge, displacement, volume
    MOMENTUM = auto() # e.g., momentum, angular momentum
    POWER = auto() # e.g., power, energy rate
    ENERGY = auto() # e.g., energy, work
    SIGNAL = auto() # A signal represents one arbitrary variable of time that may also be an effort or a flow, but not necessarily
    CONSTANT = auto() # physical constants or parameters.
@dataclass(eq=False)
class BGPhyQuantity:
    """Metadata for a bond graph variable."""
    id: str
    type: BGVariable = BGVariable.EFFORT
    physical_quantity: PhysicalQuantity | None = None

@dataclass
class Equation:
    """Represents a single equation: y = f(x) with a description."""
    dependent_symbol: str = ""
    infix_rhs: str = ""
    voi: str = "" # only for ode
    expression: str = ""  # Could be a string for now, or a sympy.Expr in a real solver
    description: str = ""

@dataclass(eq=False)
class Port:
    """Represents one typed connection point on a component."""
    label: str # the ports of each component are uniquely labelled
    component: Component = field(repr=False) # Prevents Infinite Recursion Crashing
    type: PortType = PortType.POWER_PORT
    storage: StorageType | None = None # Only relevant for storage elements; otherwise None
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
    def effort(self) -> BGPhyQuantity:
        """Returns this port's effort description, if assigned."""
        # Returns the stored symbol, or None if it hasn't been set yet.
        return BGPhyQuantity(f'e_{self.name.replace('.', '_')}',BGVariable.EFFORT, getattr(self, '_effort', None))  

    @effort.setter
    def effort(self, physical: PhysicalQuantity) -> None:
        """Sets this port's effort description."""
        self._effort = physical

    @property
    def flow(self) -> BGPhyQuantity:
        """Returns this port's flow description, if assigned."""
        return BGPhyQuantity(f'f_{self.name.replace('.', '_')}',BGVariable.FLOW, getattr(self, '_flow', None))  

    @flow.setter
    def flow(self, physical: PhysicalQuantity) -> None:
        """Sets this port's flow description."""
        self._flow = physical

    @property
    def quantity(self) -> BGPhyQuantity:
        """Returns this port's quantity description, if assigned."""
        return BGPhyQuantity(f'q_{self.name.replace('.', '_')}',BGVariable.QUANTITY, getattr(self, '_quantity', None))  

    @quantity.setter
    def quantity(self, physical: PhysicalQuantity) -> None:
        """Sets this port's quantity description."""
        self._quantity = physical
    
    @property
    def momentum(self) -> BGPhyQuantity:
        """Returns this port's momentum description, if assigned."""
        return BGPhyQuantity(f"p_{self.name.replace('.', '_')}", BGVariable.MOMENTUM, getattr(self, '_momentum', None))  
    
    @momentum.setter
    def momentum(self, physical: PhysicalQuantity) -> None:
        """Sets this port's momentum description."""
        self._momentum = physical

    @property
    def signal(self) -> BGPhyQuantity:
        """Returns this port's signal description, if assigned."""
        return BGPhyQuantity(f"s_{self.name.replace('.', '_')}", BGVariable.SIGNAL, getattr(self, '_signal', None))  
    
    @signal.setter
    def signal(self, physical: PhysicalQuantity) -> None:
        """Sets this port's signal description."""
        self._signal = physical
@dataclass(eq=False)
class Bond:
    """Connects two ports and owns their shared causality assignment."""
    source: Port
    target: Port
    type: ConnectionType = ConnectionType.POWER_BOND

    @property
    def name(self) -> str:
        """Returns the source-to-target bond identifier."""
        return f"{self.source.name}--{self.target.name}"
    
    def __post_init__(self) -> None:
        """Validates endpoints and atomically attaches the bond to both ports."""
        self.source._attach_bond(self)
        self.target._attach_bond(self)

    def disconnect(self) -> None:
        """Safely severs the bidirectional link between the bond and its ports."""
        self.source._detach_bond(self)
        self.target._detach_bond(self)

    def validate_connection(self) -> None:
        """Validates the bond's source and target ports."""
        Bond.validate(self.source, self.target, self.type)

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

    def validate_causality(self) -> bool:
        source = self.source.causality
        target = self.target.causality
    
        if source is None and target is None:
            raise ValueError(
                f"Bond '{self.name}' has unassigned causality for both endpoints."
            )
    
        if source is None or target is None:
            raise ValueError(
                f"Bond '{self.name}' has partially assigned causality."
            )
    
        if source == target:
            raise ValueError(
                f"Bond '{self.name}' has conflicting causality: "
                "both endpoints have the same causality."
            )
        return True
    
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
        
    @staticmethod
    def validate(source: Port, target: Port, connection_type: ConnectionType = ConnectionType.POWER_BOND) -> None:

        if source is target:
            raise ValueError(f"Cannot create a bond from port '{source.name}' to itself.")

        if source.component is target.component:
            raise ValueError("Cannot create a bond between ports on the same component.")

        if source.bond is not None:
            raise ValueError( f"Source port {source.name} is already connected to a bond.")

        if target.bond is not None:
            raise ValueError(f"Target port {target.name} is already connected to a bond.")

        if connection_type == ConnectionType.POWER_BOND:
           
            if source.type is not PortType.POWER_PORT:
                raise ValueError( f"Source port {source.name} is not a valid power port.")

            if target.type is not PortType.POWER_PORT:
                raise ValueError(
                    f"Target port {target.name} is not a valid power port."
                )
        elif connection_type == ConnectionType.SIGNAL_BOND:
            if source.type is not PortType.SIGNAL_PORT:
                raise ValueError( f"Source port {source.name} is not a valid signal port.")

            if target.type is not PortType.SIGNAL_PORT:
                raise ValueError( f"Target port {target.name} is not a valid signal port.")
        else:
            raise ValueError(
                f"Unsupported connection type: {connection_type!r}."
            )

@dataclass(eq=False)
class Component:
    """Models a bond-graph element, its ports, parameters, states, and equations."""
    name: str
    type: ComponentType | str = ComponentType.CUSTOM
    domain: Domain | str = Domain.ABSTRACT
    non_invertible: bool = False # If True, the component has any constitutive relationship that cannot be algebraically inverted to solve for either effort or flow.
    # Optional parameters strictly for ComponentType.CUSTOM
    num_power_ports: int = 1 # Number of power ports for ComponentType.CUSTOM, must be >= 0
    num_signal_ports: int = 0 # Number of signal ports for ComponentType.CUSTOM and Block Diagram elements, must be >= 0
    ports: dict[str, Port] = field(default_factory=dict, repr=False, init=False)
    _available_ports: deque[Port] = field(default_factory=deque, repr=False,   init=False) # track which ports are available for new bonds  
    _next_port_number: int = field(default=1, repr=False, init=False) # only used for junctions, to auto-label new ports
    bonds: set[Bond] = field(default_factory=set, repr=False, init=False) # register all bonds connected to this component, for quick lookup and deletion 
    constitutive_equations: list[str] = field(default_factory=list, repr=False, init=False) # store any constitutive equations for this component
    parameters: set[BGPhyQuantity] | None = field(default=None, repr=False, init=False) # store any parameters for this component

    def __post_init__(self) -> None:
        """Creates default ports based on the component type and requested count."""
        # 1-Port Elements
        if self.type in (ComponentType.R, ComponentType.SE, ComponentType.SF,ComponentType.C,ComponentType.I):
            storageType = StorageType.C_TYPE if self.type == ComponentType.C else StorageType.I_TYPE if self.type == ComponentType.I else None
            p1=self._add_port("p1",storage=storageType)
            self._available_ports.append(p1) # For 1-port elements, the single port is always available for bonding
        # 2-Port Elements
        elif self.type in (ComponentType.TF, ComponentType.GY,ComponentType.IC):
            storageType = StorageType.C_TYPE if self.type == ComponentType.IC else None
            p1=self._add_port("p1",storage=storageType)
            storageType = StorageType.I_TYPE if self.type == ComponentType.IC else None
            p2=self._add_port("p2",storage=storageType)
            self._available_ports.extend([p1, p2])
        elif self.type in (ComponentType.MSE, ComponentType.MSF,ComponentType.MC, ComponentType.MI):
            storageType = StorageType.C_TYPE if self.type == ComponentType.MC else StorageType.I_TYPE if self.type == ComponentType.MI else None
            p1=self._add_port("p1",storage=storageType)
            p2=self._add_port("mod", type=PortType.SIGNAL_PORT)
            self._available_ports.extend([p1, p2])
        # Reaction Elements
        elif self.type == ComponentType.Re:
            p1=self._add_port("p1",fixed_causality=True)
            p2=self._add_port("p2", fixed_causality=True)
            self._available_ports.extend([p1, p2])
            self.non_invertible = True # Reactions are generally non-invertible due to their nonlinear constitutive relationships
        # 3-Port Elements
        elif self.type in (ComponentType.MIC,ComponentType.MTF, ComponentType.MGY):
            p1=self._add_port("p1")
            p2=self._add_port("p2")
            p3=self._add_port("mod", type=PortType.SIGNAL_PORT)
            self._available_ports.extend([p1, p2, p3])
        elif self.type == ComponentType.Re_GHK:
            p1=self._add_port("p1", fixed_causality=True)
            p2=self._add_port("p2", fixed_causality=True)
            p3=self._add_port("mod", type=PortType.SIGNAL_PORT)
            self._available_ports.extend([p1, p2, p3])
            self.non_invertible = True # Modulated storage elements are generally non-invertible due to their nonlinear constitutive relationships
        elif self.type in (ComponentType.ZERO, ComponentType.ONE, ComponentType.XZERO, ComponentType.XONE):
            pass # Junctions dynamically allocate ports as needed; no default ports are created.
        elif self.type == ComponentType.BLOCK:
            for i in range(1, self.num_signal_ports + 1):
                p = self._add_port(f"s{i}", type=PortType.SIGNAL_PORT)
                self._available_ports.append(p)
        else: # type == ComponentType.CUSTOM or any other unrecognized type 
            # For custom components, create the specified number of power and signal ports
            for i in range(1, self.num_power_ports + 1):
                p = self._add_port(f"p{i}", type=PortType.POWER_PORT)
                self._available_ports.append(p)
            for i in range(1, self.num_signal_ports + 1):
                p = self._add_port(f"s{i}", type=PortType.SIGNAL_PORT)
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

    def add_parameter(self, name: str, physical_quantity: PhysicalQuantity ) -> None:
        """Adds a parameter to the component's parameter set."""
        if self.parameters is None:
            self.parameters = set()
        self.parameters.add(BGPhyQuantity(id=name, type=BGVariable.CONSTANT, physical_quantity=physical_quantity))

    def add_constitutive_equation(self, equation: str) -> None:
        """Adds a constitutive equation to the component's list."""
        self.constitutive_equations.append(equation)

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
        if self.type in JUNCTIONS:
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
                f"Component '{self.name}' of type '{self.type}' does not support dynamic port allocation."
            )   

    def get_available_ports(self) -> list[Port]:
        """Returns a list of currently unconnected ports."""
        # This method is only relevant for non-junction components. For junctions, use `allocate_port()` to get a free port or create a new one.
        return list(self._available_ports)

    def set_fixed_causality(self, port_label: str, causality_value: bool) -> None:
        """
        Assigns a fixed causality to a specific port after component creation.
        
        Args:
            port_label: The label of the port in the self.ports dictionary.
            causality_value: The causality state to assign .
        """
        if not self.ports:
            raise ValueError(f"Component '{self.name}' has no ports initialized.")
            
        try:
            self.ports[port_label].fixed_causality = causality_value
        except KeyError:
            raise KeyError(
                f"Port label '{port_label}' not found in component '{self.name}'."
            )
        self.non_invertible = True # Mark the component as non-invertible if a fixed causality is set on any port
class BondGraph:
    """Owns a connected set of components"""

    def __init__(self, name: str = "bond_graph") -> None:
        """Initializes an empty named graph."""
        self.name = name
        self.components: dict[str, Component] = {} # Mapping of component names to Component objects
        # Insertion-ordered mapping gives O(1) bond membership/deletion
        # while retaining deterministic iteration order.
        self._bonds: dict[Bond, None] = {}
        self.physical_constants: set[BGPhyQuantity] | None = None # Global physical constants for the bond graph
        self.equations: list[Equation] | None = None # Store any equations for the bond graph
     
    @property
    def bonds(self):
        """Insertion-ordered, set-like view of all bonds."""
        return self._bonds.keys()

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
                if arg.type in JUNCTIONS:
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

    def set_fixed_causality(self, comp_arg: Component | str, port_label: str, causality_value: bool) -> None:
        """Sets a fixed causality for a specific port on a component."""
        comp_name = comp_arg.name if isinstance(comp_arg, Component) else comp_arg
        comp = self.components.get(comp_name)
        if not comp:
            warnings.warn(f"Component '{comp_name}' not found in bond graph.")
            return
        try:
            comp.set_fixed_causality(port_label, causality_value)
        except KeyError as e:
            warnings.warn(str(e))

    def set_domain(self, arg: Component| Port| str, domain_value: Domain | str) -> None:
        """Sets the domain for a specific component."""
        if isinstance(arg, Component):
            arg.domain = domain_value
        elif isinstance(arg, Port):
            arg.domain = domain_value
        elif isinstance(arg, str):
            resolved = self._resolve_string(arg)
            if isinstance(resolved, Component) or isinstance(resolved, Port):
                resolved.domain = domain_value
            else:
                warnings.warn(f"Could not resolve '{arg}' to a valid component or port for domain assignment.")
        else:
            warnings.warn(f"Invalid argument type: {type(arg)}. Expected Component, Port, or str.")

    def add_physical_constant(self, name: str, physical_quantity: PhysicalQuantity) -> None:
        """Adds a global physical constant to the bond graph."""
        if self.physical_constants is None:
            self.physical_constants = set()
        self.physical_constants.add(BGPhyQuantity(id=name, type=BGVariable.CONSTANT, physical_quantity=physical_quantity))

class DomainRefiner:
    """Decorates abstract Bond Graph components with physical domain knowledge."""
    
    def __init__(self, catalog_path: str) -> None:
        with open(catalog_path, 'r') as f:
            self.catalog = json.load(f)

    def refine_component(self, component: Component, domain: Domain, template_id: str | None = None) -> None:
        """Applies domain variables, parameters, and equations to an existing component."""
        component.domain = domain
        domain_data = self.catalog.get(domain.name, {})

        if not template_id:
            template_id = getattr(component.type, 'name', str(component.type))
            
        comp_metadata = domain_data.get("components", {}).get(template_id, {})
        port_domain_overrides = comp_metadata.get("port_domains", {})
        multiports = (len(component.ports) > 1)      
        # 1. Map domain variables to ports (handling multi-domain overrides)
        for port_label, port in component.ports.items():
            # Check if this specific port has a designated domain in the JSON
            override_domain_str = port_domain_overrides.get(port_label)
            
            if override_domain_str:
                # Resolve the string to your Domain enum
                port_domain_enum = getattr(Domain, override_domain_str, Domain.ABSTRACT)
                port.domain = port_domain_enum
                
                # Fetch the correct variable definitions from the overarching catalog
                port_domain_data = self.catalog.get(override_domain_str, {})
                domain_vars = port_domain_data.get("domain_variables", {})
            else:
                # Fallback to the component's primary domain
                port.domain = domain
                domain_vars = domain_data.get("domain_variables", {})
            # Apply the variables via the @property setters
            for var_key in ["effort", "flow", "quantity", "momentum", "signal"]:
                if var_key in domain_vars:
                    # Copy before mutating so repeated calls don't keep appending to the shared catalog entry.
                    var_data = dict(domain_vars[var_key])
                    if "symbol" in var_data:
                        suffix = port.name.replace('.', '_') if multiports else component.name
                        var_data["symbol"] = f"{var_data['symbol']}_{suffix}"
                    pq = _pq_from_serializable(var_data)
                    if pq is not None:
                        setattr(port, var_key, pq)
        # 2. Apply Equations & Parameters
        if comp_metadata:
            for eq in comp_metadata.get("constitutive_equations", []):
                if eq not in component.constitutive_equations:
                    component.add_constitutive_equation(eq)
            
            for p_name, p_data in comp_metadata.get("parameters", {}).items():
                # Copy before mutating so repeated calls don't keep appending to the shared catalog entry.
                p_data = dict(p_data)
                if "symbol" in p_data:
                    p_data["symbol"] = f"{p_data['symbol']}_{component.name}"
                pq = _pq_from_serializable(p_data)
                if pq is not None:
                    component.add_parameter(p_name, pq)
                else:
                    warnings.warn(f"Invalid physical quantity data for parameter '{p_name}' in component '{component.name}'.")

    def refine_graph(self, bg: BondGraph, refinement_map: dict[str, tuple[Domain, str]]) -> None:
        """Batch refines a whole graph and extracts global parameters."""
    
        # 1. Apply global parameters to the graph first
        for domain, _ in refinement_map.values():
            domain_data = self.catalog.get(domain.name, {})
            for g_name, g_data in domain_data.get("physical_constants", {}).items():
                if g_name not in bg.physical_constants:                    
                    g_pq = _pq_from_serializable(g_data)
                    if g_pq is not None:
                        bg.add_physical_constant(g_name, g_pq)
                    else:
                        warnings.warn(f"Invalid physical quantity data for global constant '{g_name}' in domain '{domain.name}'.")
                else:
                    warnings.warn(f"Global constant '{g_name}' already exists in the bond graph; skipping duplicate.")
        # 2. Refine individual components
        for comp_name, (domain, template_id) in refinement_map.items():
            comp = bg.components.get(comp_name)
            if comp:
                self.refine_component(comp, domain, template_id)


def _pq_to_serializable(physical_quantity: PhysicalQuantity) -> dict:
    """Converts a PhysicalQuantity dataclass to a dict, passing through other values unchanged."""
    return asdict(physical_quantity) if is_dataclass(physical_quantity) else physical_quantity

def _pq_from_serializable(data: dict) -> PhysicalQuantity | None:
    """Reconstructs a PhysicalQuantity from an imported dict. 'value' is optional; the rest are required."""
    if isinstance(data, dict) and {"description", "symbol", "units"}.issubset(data.keys()):
        return PhysicalQuantity(**data)
    return None

def _eq_to_serializable(equation: Equation) -> dict:
    """Converts an Equation dataclass to a dict, passing through other values unchanged."""
    return asdict(equation) if is_dataclass(equation) else equation

def _eq_from_serializable(data: dict) -> Equation | None:
    """Reconstructs an Equation from an imported dict. 'expression' is required; the rest are optional."""
    if isinstance(data, dict) and {"dependent_symbol", "infix_rhs", "voi", "expression"} .issubset(data.keys()):
        return Equation(**data)
    return None

def exportBG(bg: BondGraph,json_file: str) -> None:
    """Serializes a BondGraph object and all state variables to a JSON string."""
    data = {
        "name": bg.name,
        "components": [],
        "bonds": []
    }
    if bg.equations is not None:
        data["equations"] = [_eq_to_serializable(eq) for eq in bg.equations]

    if bg.physical_constants is not None:
        data["physical_constants"] = {}
        for pq in bg.physical_constants:
            data["physical_constants"][pq.id] = _pq_to_serializable(pq.physical_quantity) if pq.physical_quantity is not None else None
    # 1. Export Components and Ports
    for comp in bg.components.values():
        c_data = {
            "name": comp.name,
            "type": comp.type.name if isinstance(comp.type, ComponentType) else comp.type,
            "domain": comp.domain.name if isinstance(comp.domain, Domain) else comp.domain,
            "non_invertible": comp.non_invertible,
            "num_power_ports": comp.num_power_ports,
            "num_signal_ports": comp.num_signal_ports,
            "ports": []
        }
        if len(comp.constitutive_equations) > 0:
            c_data["constitutive_equations"] = comp.constitutive_equations
        if comp.parameters is not None:
            c_data["parameters"] = {}
            for  pq in comp.parameters:
                c_data["parameters"][pq.id] = _pq_to_serializable(pq.physical_quantity) if pq.physical_quantity is not None else None

        for port_label, port in comp.ports.items():
            p_data = {
                "label": port_label,
                "type": port.type.name,
                "storage": port.storage.name if port.storage else None,
                "domain": port.domain.name if isinstance(port.domain, Domain) else port.domain,
                "fixed_causality": port.fixed_causality,
                "causality": port.causality
                }
            if port.effort.physical_quantity is not None:
                p_data["effort"] = _pq_to_serializable(port.effort.physical_quantity)
            if port.flow.physical_quantity is not None:
                p_data["flow"]= _pq_to_serializable(port.flow.physical_quantity)
            if port.quantity.physical_quantity is not None:
                p_data["quantity"] = _pq_to_serializable(port.quantity.physical_quantity)
            if port.momentum.physical_quantity is not None:
                p_data["momentum"] =  _pq_to_serializable(port.momentum.physical_quantity)
            if port.signal.physical_quantity is not None:
                p_data["signal"]=_pq_to_serializable(port.signal.physical_quantity)

            c_data["ports"].append(p_data)
            
        data["components"].append(c_data)

    # 2. Export Bonds
    for bond in bg.bonds:
        data["bonds"].append({
            "source": bond.source.name,  # Uses Port.name (Component.label)
            "target": bond.target.name,  # Uses Port.name (Component.label)
            "type": bond.type.name
        })

    with open(json_file, "w") as f:
        json.dump(data, f, indent=4)

def importBG(json_file: str) -> BondGraph:
    """Reconstructs a BondGraph object from a JSON file."""
    with open(json_file, "r") as f:
        data = json.load(f)
    bg = BondGraph(name=data.get("name", "Imported_BG"))
    if "equations" in data:
        bg.equations = []
        for eq in data["equations"]:
            eq_obj = _eq_from_serializable(eq)
            if eq_obj is not None:
                bg.equations.append(eq_obj)
    if "physical_constants" in data:
        bg.physical_constants = set()
        for name, pq in data["physical_constants"].items():
            pq_obj = _pq_from_serializable(pq)
            if pq_obj is not None:
                bg.add_physical_constant(name, pq_obj)
    # 1. Reconstruct Components and Ports
    for c_data in data.get("components", []):
        c_type_val = c_data["type"]
        domain_val = c_data["domain"]
        
        c_type = getattr(ComponentType, c_type_val, c_type_val)
        domain = getattr(Domain, domain_val, domain_val)
        
        comp = bg.add_component(
            c_data["name"], 
            type=c_type, 
            domain=domain,
            num_power_ports=c_data.get("num_power_ports", 1),
            num_signal_ports=c_data.get("num_signal_ports", 0)
        )
        comp.non_invertible = c_data.get("non_invertible", False)
        if "constitutive_equations" in c_data:
            comp.constitutive_equations = c_data["constitutive_equations"]
        if "parameters" in c_data:
            comp.parameters = set()
            for name, pq in c_data["parameters"].items():
                pq_obj = _pq_from_serializable(pq)
                if pq_obj is not None:
                    comp.add_parameter(name, pq_obj)
        # Restore precise port states (crucial for junctions which do not auto-generate ports)
        for p_data in c_data.get("ports", []):
            label = p_data["label"]
            p_domain = getattr(Domain, p_data.get("domain")) if p_data.get("domain") in Domain.__members__ else p_data.get("domain", Domain.ABSTRACT)

            # If __post_init__ didn't create the port (e.g., Junctions), create it now
            if label not in comp.ports:
                p_type = getattr(PortType, p_data.get("type", "POWER_PORT"))
                comp._add_port(label, type=p_type, domain=p_domain)
            # Keep the _next_port_number updated for junctions
            if label.startswith("p"):
                try:
                    port_num = int(label[1:])
                    if port_num >= comp._next_port_number:
                        comp._next_port_number = port_num + 1
                except ValueError:
                    pass                   
            # Apply state variables
            port = comp.ports[label]
            port.domain = p_domain  # Ports created eagerly by __post_init__ also need the restored domain
            port.fixed_causality = p_data.get("fixed_causality")
            port.causality = p_data.get("causality")

            # Get the physical quantities for each BGVariable type
            for var_type in ["effort", "flow", "quantity", "momentum", "signal"]:
                var_data = p_data.get(var_type)
                if var_data:
                    # port.effort/flow/... are properties backed by _effort/_flow/...;
                    # assign through the setter, not the throwaway getter object.
                    if _pq_from_serializable(var_data) is not None:
                        setattr(port, var_type, _pq_from_serializable(var_data))
            
            # Re-register port availability
            if not port.bond and port not in comp._available_ports:
                comp._available_ports.append(port)

    # 2. Reconstruct Bonds
    for b_data in data.get("bonds", []):
        b_type = getattr(ConnectionType, b_data.get("type", "POWER_BOND"))
        
        # Uses _resolve_string natively to map "Component.label" to the correct port
        bg.add_bond(
            source=b_data["source"],
            target=b_data["target"],
            type=b_type
        )

    return bg

def importBG_pq(bg: BondGraph, json_file: str) -> BondGraph:
    """Imports only the physical quantities of a BondGraph from a JSON file, including constants,
    component parameters, and port variables."""
    with open(json_file, "r") as f:
        data = json.load(f)

    if "physical_constants" in data:
        bg.physical_constants = set()
        for name, pq in data["physical_constants"].items():
            pq_obj = _pq_from_serializable(pq)
            if pq_obj is not None:
                bg.add_physical_constant(name, pq_obj)

    for c_data in data.get("components", []):
        comp = bg.components.get(c_data["name"])
        if not comp:
            warnings.warn(f"Component '{c_data['name']}' not found in bond graph; skipping parameter and port PQ import.")
            continue
        if "parameters" in c_data:
            comp.parameters = set()
            for name, pq in c_data["parameters"].items():
                pq_obj = _pq_from_serializable(pq)
                if pq_obj is not None:
                    comp.add_parameter(name, pq_obj)
        for p_data in c_data.get("ports", []):
            label = p_data["label"]
            if label not in comp.ports:
                warnings.warn(f"Port '{label}' not found on component '{comp.name}'; skipping PQ import for this port.")
                continue
            else:
                port = comp.ports[label]
            for var_type in ["effort", "flow", "quantity", "momentum", "signal"]:
                var_data = p_data.get(var_type)
                if var_data:
                    pq_obj = _pq_from_serializable(var_data)
                    if pq_obj is not None:
                        setattr(port, var_type, pq_obj)
    return bg

if __name__ == "__main__": 
    # Construct the system: Mass (I), Spring (C), Damper (R), Force Source (SE)
    bg = BondGraph("Mass_Spring_Damper")
    # Add components
    se = bg.add_component("SE_Force", type=ComponentType.SE)
    j1 = bg.add_component("J1", type=ComponentType.ONE)
    mass = bg.add_component("I_Mass", type=ComponentType.I)
    spring = bg.add_component("C_Spring", type=ComponentType.C)
    damper = bg.add_component("R_Damper", type=ComponentType.R)
    # Connect components via bonds
    bg.add_bond(se, j1)
    bg.add_bond(j1, mass)
    bg.add_bond(j1, spring)
    bg.add_bond(j1, damper)

    # 1. Refine the graph with domain knowledge
    refiner = DomainRefiner("domain_catalog.json")
    refiner.refine_graph(bg, {
        "I_Mass": (Domain.MECHANICAL_TRANSLATIONAL, "I"),
        "C_Spring": (Domain.MECHANICAL_TRANSLATIONAL, "C"),
        "R_Damper": (Domain.MECHANICAL_TRANSLATIONAL, "R"),
        "SE_Force": (Domain.MECHANICAL_TRANSLATIONAL, "Se")
    })

    # 1. Export unassigned or partially assigned graph
    exportBG(bg, "mass_spring_damper.json")
    restored_bg = importBG("mass_spring_damper.json")
    exportBG(restored_bg, "restored_mass_spring_damper.json")
    imported_bg = importBG_pq(restored_bg, "mass_spring_damper.json")
    exportBG(imported_bg, "imported_mass_spring_damper_pq.json")
