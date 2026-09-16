import xml.etree.ElementTree as ET
from utilities import infix_to_mathml
from xml.dom import minidom
from defineBG import BondGraph,PhysicalQuantity, importBG

CellMLV1_namespaces = {
        'cellml': "http://www.cellml.org/cellml/1.1#",  # CellML namespace
        'xlink': "http://www.w3.org/1999/xlink",  # XLink namespace
        'math': "http://www.w3.org/1998/Math/MathML"  # MathML namespace
    }
defUnit=["ampere","becquerel","candela","celsius","coulomb","dimensionless","farad","gram","gray","henry",
    "hertz","joule","katal","kelvin","kilogram","liter","litre","lumen","lux","meter","metre","mole",
    "newton","ohm","pascal","radian","second","siemens","sievert","steradian","tesla","volt","watt","weber"]

# modified from https://github.com/CellDL/bondgraph-tools/blob/main/bondgraph/bondgraph/cellml/__init__.py

def cellml_element(tag: str, *args, **attributes) -> ET.Element:
    return ET.Element(tag, *args, **attributes)

def cellml_subelement(parent: ET.Element, tag: str, *args, **attributes) -> ET.Element:
    return ET.SubElement(parent, tag, *args, **attributes)
class CellMLVariable:    
    def __init__(self, name: str, units: str):
        self.__name = name
        self.__units = units
        self._public_interface = None
        self._private_interface = None
        self.__initial_value = None

    def set_initial_value(self, value: str | None):
        self.__initial_value = value
    def set_public_interface(self, public_interface: str | None):
        self._public_interface = public_interface
    def set_private_interface(self, private_interface: str | None):
        self._private_interface = private_interface

    def get_element(self) -> ET.Element:
        element = cellml_element('variable', name=self.__name, units=self.__units)
        if self.__initial_value is not None:
            element.attrib['initial_value'] = f'{self.__initial_value}'
        if self._public_interface is not None:
            element.attrib['public_interface'] = self._public_interface
        if self._private_interface is not None:
            element.attrib['private_interface'] = self._private_interface
        return element

    @property
    def name(self):
    #==============
        return self.__name
    @property
    def units(self):
    #==============
        return self.__units
    @property
    def initial_value(self):
    #==============
        return self.__initial_value
    @property
    def public_interface(self):
    #==============
        return self._public_interface
    @property
    def private_interface(self):
    #==============
        return self._private_interface


def create_cellmlV1_rootET(model_name, cellml_prefix=True):
    """Create a CellML V1.x model as an xml.etree.ElementTree

    Parameters
    ----------
    model_name : str
        The name of the model
    cellml_prefix : bool, optional
        Whether to include the 'cellml' prefix in the xmlns attribute, by default True

    Returns
    -------
    xml.etree.ElementTree with the root model element

    """
    # Create the root model element without auto-generated prefixes
    if cellml_prefix:
        model_attrs = {
            'name': model_name,
            'xmlns': CellMLV1_namespaces['cellml'],  
            'xmlns:cellml': CellMLV1_namespaces['cellml'],  # Explicitly add the prefix for CellML
            'xmlns:xlink': CellMLV1_namespaces['xlink']  # Explicitly add the prefix for XLink
        }
    else: # later on, we will add the cellml prefix to the math element
        model_attrs = {
            'name': model_name,
            'xmlns': CellMLV1_namespaces['cellml'],  
            'xmlns:xlink': CellMLV1_namespaces['xlink']  # Explicitly add the prefix for XLink
        }
    
    return cellml_element('model', **model_attrs)

def write_cellmlV1 (model,model_file):
    """Write the model to a CellML V1.x file
    
    Parameters
    ----------
    model : xml.etree.ElementTree
        The model to be written to a file
    model_file : str
        The file path of the CellML file

    Returns
    -------
    None      
    """
    # Convert to string using ElementTree, then parse with minidom for pretty printing
    # Ensure model is an Element, not an ElementTree
    if isinstance(model, ET.ElementTree):
        model = model.getroot()
    rough_string = ET.tostring(model, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    pretty_string = reparsed.toprettyxml(indent="    ")

    # Write the pretty-printed XML to file
    with open(model_file, 'w') as output_file:
        output_file.write(pretty_string)    

def read_cellmlV1 (model_file):
    """Read the model from a CellML V1.x file
    
    Parameters
    ----------
    model_file : str
        The file path of the CellML file

    Returns
    -------
    model : xml.etree.ElementTree
        The model read from the file
    """

    # Use the namespace dictionary to parse the CellML file
    ET.register_namespace('cellml', CellMLV1_namespaces['cellml'])
    ET.register_namespace('xlink', CellMLV1_namespaces['xlink'])
    ET.register_namespace('math', CellMLV1_namespaces['math'])    
   # ET.register_namespace('', CellMLV1_namespaces['cellml'])
    model = ET.parse(model_file).getroot()
    # Remove the 'cellml' namespace prefix from element tags
    def remove_namespace_prefix(element, prefix):
        if element.tag.startswith(f'{{{CellMLV1_namespaces[prefix]}}}'):
            element.tag = element.tag.replace(f'{{{CellMLV1_namespaces[prefix]}}}', '')
        for child in element:
            remove_namespace_prefix(child, prefix)

    # Update the default namespace in the model and remove the 'cellml' prefix
    model.attrib['xmlns'] = CellMLV1_namespaces['cellml']
    remove_namespace_prefix(model, 'cellml')  
    return model

def units_import(model_ET, units_Set,units_file):
    """Add units import to a CellML V1.x model as ET. ElementTree

    Parameters
    ----------
    model_ET : xml.etree.ElementTree
        The model to which the units import will be added
    units_Set : set
        A set of units to be added
    units_file : str
        The file path of the units file

    Returns
    -------
    None

    Side effects
    ------------
    The model_ET is modified in place

    """

    if len(units_Set) == 0:
        return
    units_import=cellml_subelement(model_ET, 'import', {'xlink:href': units_file})
    for units_name in units_Set:
        if units_name not in defUnit:
            units_import_i = cellml_subelement(units_import, 'units', {'name': units_name, 'units_ref': units_name})


def pq2CellMLVariable(pq: PhysicalQuantity) -> CellMLVariable:
    """Convert a PhysicalQuantity to a CellMLVariable."""
    cellml_var = CellMLVariable(pq.symbol, pq.units)
    if pq.value is not None:
        cellml_var.set_initial_value(str(pq.value))    
    return cellml_var

def BG2CellMLV1(bg: BondGraph):

    module_name = bg.name
    params = []
    vars = []
    units_set=set()
    if bg.physical_constants is not None:
        for const in bg.physical_constants:
            if const.physical_quantity is not None:
                cellml_var = pq2CellMLVariable(const.physical_quantity)
                params.append(cellml_var.get_element())
                units_set.add(const.physical_quantity.units)
    for comp in bg.components.values():
        if comp.parameters is not None:
            for param in comp.parameters:
                if param.physical_quantity is not None:
                    cellml_var = pq2CellMLVariable(param.physical_quantity)
                    params.append(cellml_var.get_element())
                    units_set.add(param.physical_quantity.units)
        for port in comp.ports.values():
            if port.effort.physical_quantity is not None:
                cellml_var = pq2CellMLVariable(port.effort.physical_quantity)
                vars.append(cellml_var.get_element())
                units_set.add(port.effort.physical_quantity.units)
            if port.flow.physical_quantity is not None:
                cellml_var = pq2CellMLVariable(port.flow.physical_quantity)
                vars.append(cellml_var.get_element())
                units_set.add(port.flow.physical_quantity.units)
            if port.quantity.physical_quantity is not None:
                cellml_var = pq2CellMLVariable(port.quantity.physical_quantity)
                vars.append(cellml_var.get_element())
                units_set.add(port.quantity.physical_quantity.units)
            if port.momentum.physical_quantity is not None:
                cellml_var = pq2CellMLVariable(port.momentum.physical_quantity)
                vars.append(cellml_var.get_element())
                units_set.add(port.momentum.physical_quantity.units)
            if port.signal.physical_quantity is not None:
                cellml_var = pq2CellMLVariable(port.signal.physical_quantity)
                vars.append(cellml_var.get_element())
                units_set.add(port.signal.physical_quantity.units)


    module_ET=create_cellmlV1_rootET(module_name,cellml_prefix=False)
    # register the CellML namespace    
    ET.register_namespace('cellml', CellMLV1_namespaces['cellml'])
    # Create the MathML element with the correct namespace
    mathml_element = cellml_element('math', {'xmlns': CellMLV1_namespaces['math']})   
    module_component = cellml_subelement(module_ET, 'component', {'name': module_name})
    if bg.equations is not None:
        for eq in bg.equations:
            mmathml_string = infix_to_mathml(eq.dependent_symbol, eq.infix_rhs, eq.voi)
            # Parse the MathML string
            try:
                # Convert the generated math string into an XML element
                math_content = ET.fromstring(mmathml_string)
                if len(math_content) > 0:
                    mathml_element.append(math_content[0])  # Append the first child, which is the actual content
            except ET.ParseError as e:
                print(f"Error parsing MathML: {e}")
            if eq.voi:
                # Add the variable of integration to the module component
                if eq.voi not in [var.attrib['name'] for var in vars]:
                    voi_var = CellMLVariable(eq.voi, "second")  # Assuming VOI is in seconds
                    vars.append(voi_var.get_element())
   
    # Add parameters to the module component
    for param_var in params:
        module_component.append(param_var)
    # Add variables to the module component
    for var in vars:
        module_component.append(var)  
    # Append the <math> element to the component
    module_component.append(mathml_element)
    # Add the units import to the module
    units_set.difference_update(defUnit)
    units_import(module_ET, units_set, 'units.cellml')

    return  module_ET
                     

if __name__ == "__main__": 
    
   bg = importBG("mass_spring_damper_equations.json")
   model_ET=BG2CellMLV1(bg)
   write_cellmlV1(model_ET, "mass_spring_damper.cellml")