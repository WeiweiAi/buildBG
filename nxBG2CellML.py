from webbrowser import get
from build_nxBG import getPowerPorts, load_nxBG_json
import xml.etree.ElementTree as ET
from utilities import infix_to_mathml
import copy
from xml.dom import minidom
from pathlib import Path
import libcellml


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

def import_comp(model_ET, comp_name_new, comp_name, model_file, imported_model):
    """Add a component import to a CellML V1.x model as ET. ElementTree

    Parameters
    ----------
    model_ET : xml.etree.ElementTree or xml.etree.Element
        The model to which the import will be added
    comp_name_new : str
        The component name in the new model
    comp_name : str
        The component name of the imported model
    model_file : str
        The file path of the imported model
    imported_model : xml.etree.ElementTree or xml.etree.Element
        The imported model

    Returns
    -------
    None

    Side effects
    ------------
    The model_ET is modified in place

    """
    
    import_component= None
    # get the component element from the imported model
    for component in imported_model.findall('component'):
        if component.attrib['name']==comp_name:
            component_element=copy.deepcopy(component)
            model_import = cellml_subelement(model_ET, 'import', {'xlink:href': model_file})
            import_component=cellml_subelement(model_import, 'component', {'name': comp_name_new, 'component_ref': comp_name})
            return component_element  
    if import_component is None:
        print('The component is not found in the imported model')
        return None

def map_components(model, component1, component2):
    """Map the variables of two components in a CellML V1.x model as ET. ElementTree

       The mapping is based on the name of the variables and also the units of the variables
       Component2 dominates the public_interface of the variables if the public_interface is not defined in component1

    Parameters
    ----------
    model : xml.etree.ElementTree or xml.etree.Element
        The model to which the mapping will be added
    component1 : xml.etree.ElementTree or xml.etree.Element
        The first component to be mapped
    component2 : xml.etree.ElementTree or xml.etree.Element
        The second component to be mapped
    Returns
    -------
    None

    Side effects
    ------------
    The model_ET is modified in place
    """
    flag_component_pair=False
    component_name1=component1.attrib['name']
    component_name2=component2.attrib['name']
    component_pair={'component_1': component_name1, 'component_2': component_name2}
    for variable_1 in component1.iter('variable'):
        for variable_2 in component2.iter('variable'):
            if variable_1.attrib['name']==variable_2.attrib['name'] and variable_1.attrib['units']==variable_2.attrib['units']:
                attribute_map_variable={'variable_1': variable_1.attrib['name'], 'variable_2': variable_2.attrib['name']}
                flag_variable_pair=False
                if 'public_interface'in variable_2.attrib:
                    if 'public_interface' not in variable_1.attrib:
                        if variable_2.attrib['public_interface']=='in':
                            variable_1.attrib['public_interface']='out'
                        if variable_2.attrib['public_interface']=='out':
                            variable_1.attrib['public_interface']='in'
                        if not flag_component_pair:
                            connection=ET.SubElement(model, 'connection')
                            map_components=ET.SubElement(connection, 'map_components', component_pair)
                            flag_component_pair=True
                        map_variables=ET.SubElement(connection, 'map_variables',attribute_map_variable)
                        flag_variable_pair=True
                    else:
                        if variable_1.attrib['public_interface']!=variable_2.attrib['public_interface']:
                            if flag_component_pair:
                                for attribute_map_variables in connection.findall('map_variables'):
                                    if attribute_map_variables==attribute_map_variable:
                                        flag_variable_pair=True
                                if not flag_variable_pair:
                                    map_variables=ET.SubElement(connection, 'map_variables', attribute_map_variable)
                                    flag_variable_pair=True
                            else:
                                connection=ET.SubElement(model, 'connection')
                                map_components=ET.SubElement(connection, 'map_components', component_pair)
                                flag_component_pair=True                                    
                                map_variables=ET.SubElement(connection, 'map_variables', attribute_map_variable)
                                flag_variable_pair=True
                        else:
                            print('The public_interface of {} is the same in both components'.format(variable_2.attrib['name']))
                            pass
                else:
                    print('The public_interface of {} is not defined'.format(variable_2.attrib['name']))

def nxBG2CellMLV1(G, implementation_dict):

    """
    implemention_dict={'language':'CellML','version':1.1, 'filepath':'./data/',
                    'module_name':'SLC5_BG','module_file':'SLC5_BG.cellml',
                    'params_name':'SLC5_BG_param','params_file':'SLC5_BG_param.cellml',
                    'model_name':'SLC5_BG_run','model_file':'SLC5_BG_run.cellml',
                    'units_file':'SLC5_BG_units.cellml', 'voi':{'propertyName': 't', 'units': 'second'},
                    'observed_vars': {var_name_1, var_name_2, ...}}
    """
    units_file=implementation_dict['units_file']
    module_name=implementation_dict['module_name']
    params_name=implementation_dict['params_name']
    model_name=implementation_dict['model_name']
    module_file=implementation_dict['module_file']
    params_file=implementation_dict['params_file']
    main_file=implementation_dict['model_file']
    filepath=implementation_dict['filepath']
    observed_vars=implementation_dict['observed_vars']
    module_ET=create_cellmlV1_rootET(module_name,cellml_prefix=False)
    params_ET=create_cellmlV1_rootET(params_name) 
    model_ET=create_cellmlV1_rootET(model_name)
    # register the CellML namespace    
    ET.register_namespace('cellml', CellMLV1_namespaces['cellml'])
    # Create the MathML element with the correct namespace
    mathml_element = cellml_element('math', {'xmlns': CellMLV1_namespaces['math']})   
    module_component = cellml_subelement(module_ET, 'component', {'name': module_name})
    model_component = cellml_subelement(model_ET, 'component', {'name': model_name})
    params_component = cellml_subelement(params_ET, 'component', {'name': params_name})

    params_set=set()
    param_units_set=set()
    var_units_set=set()
    var_set=set() 
    params = []
    param_vars = []
    vars = []
    for node in G.nodes:
        if G.nodes[node]['a']=='BondElement':
            if 'modelParameter' in G.nodes[node].keys():
                for param in G.nodes[node]['modelParameter']:                    
                    if G.nodes[node]['modelParameter'][param]['propertyName'] not in params_set: 
                        var_name = G.nodes[node]['modelParameter'][param]['propertyName']
                        var_units = G.nodes[node]['modelParameter'][param]['units']
                        var_value = G.nodes[node]['modelParameter'][param]['value']
                        # put the parameters in the params component
                        param_units_set.add(var_units)
                        params_set.add(var_name)
                        cellml_var = CellMLVariable(var_name, var_units)
                        cellml_var.set_initial_value(str(var_value)) 
                        cellml_var.set_public_interface('out')                                        
                        params.append(cellml_var.get_element())
                        # put the parameters in the module component
                        var_units_set.add(var_units)
                        var_set.add(var_name)
                        cellml_var = CellMLVariable(var_name, var_units)
                        cellml_var.set_public_interface('in')
                        param_vars.append(cellml_var.get_element())
            if 'modelState' in G.nodes[node].keys():
                for state in G.nodes[node]['modelState']:
                    var_name = G.nodes[node]['modelState'][state]['propertyName']
                    var_units = G.nodes[node]['modelState'][state]['units']
                    var_init_name = G.nodes[node]['modelState'][state]['propertyName']+'_init'
                    # put the initial value as a parameter in the params component
                    param_units_set.add(var_units)
                    cellml_var = CellMLVariable(var_init_name, var_units)
                    if 'value' in G.nodes[node]['modelState'][state].keys():
                        cellml_var.set_initial_value(str(G.nodes[node]['modelState'][state]['value']))
                    else:
                        cellml_var.set_initial_value('1') # default value
                    cellml_var.set_public_interface('out')
                    params.append(cellml_var.get_element())
                    # put the state variable in the module component
                    var_units_set.add(var_units)
                    var_set.add(var_name)
                    cellml_var = CellMLVariable(var_name, var_units)
                    cellml_var.set_public_interface('in')
                    cellml_var.set_initial_value(var_init_name)
                    vars.append(cellml_var.get_element())
            if 'constitutive_eqs' in G.nodes[node].keys():
                for LHS_name, equation in G.nodes[node]['constitutive_eqs'].items():
                    mmathml_string = infix_to_mathml(LHS_name, equation[0], equation[1])
                    # Parse the MathML string
                    try:
                        # Convert the generated math string into an XML element
                        math_content = ET.fromstring(mmathml_string)
                        if len(math_content) > 0:
                            mathml_element.append(math_content[0])  # Append the first child, which is the actual content
                    except ET.ParseError as e:
                        print(f"Error parsing MathML: {e}")
            powerPorts=getPowerPorts(G,node)
            for powerPort in powerPorts:
                causality=G.nodes[powerPort]['causality']
                var_name= G.nodes[powerPort][causality]['propertyName']
                var_units= G.nodes[powerPort][causality]['units']
                cellml_var = CellMLVariable(var_name, var_units)                
                var_units_set.add(var_units)
                if cellml_var.name not in var_set:
                    vars.append(cellml_var.get_element())
                    var_set.add(var_name) 
        if G.nodes[node]['a']=='JunctionStructure' and (G.nodes[node]['subClass']=='TF' or G.nodes[node]['subClass']=='GY'):
            if 'modelParameter' in G.nodes[node].keys() and len(G.nodes[node]['modelParameter'])>1:
                for param in G.nodes[node]['modelParameter']:                    
                    if G.nodes[node]['modelParameter'][param]['propertyName'] not in params_set: 
                        var_name = G.nodes[node]['modelParameter'][param]['propertyName']
                        var_units = G.nodes[node]['modelParameter'][param]['units']
                        var_value = G.nodes[node]['modelParameter'][param]['value']
                        # put the parameters in the params component
                        param_units_set.add(var_units)
                        params_set.add(var_name)
                        cellml_var = CellMLVariable(var_name, var_units)
                        cellml_var.set_initial_value(str(var_value)) 
                        cellml_var.set_public_interface('out')                                        
                        params.append(cellml_var.get_element())
                        # put the parameters in the module component
                        var_units_set.add(var_units)
                        var_set.add(var_name)
                        cellml_var = CellMLVariable(var_name, var_units)
                        cellml_var.set_public_interface('in')
                        param_vars.append(cellml_var.get_element())

    if 'voi' in implementation_dict.keys():
        var_name = implementation_dict['voi']['propertyName']
        var_units = implementation_dict['voi']['units']
        cellml_var = CellMLVariable(var_name, var_units)
        var_units_set.add(var_units)
        vars.append(cellml_var.get_element())
    # Add parameters to the parameters component
    for param in params:
        params_component.append(param)
    units_import(params_ET, param_units_set,units_file)
    # Add parameters to the module component
    for param_var in param_vars:
        module_component.append(param_var)
    # Add variables to the module component
    for var in vars:
        module_component.append(var)  
    # Append the <math> element to the component
    module_component.append(mathml_element)
    # Add the units import to the module
    units_import(module_ET, var_units_set | param_units_set, units_file)
     
    for observable in observed_vars:
        variable_attributes = {}
        # Find the variable element with the name matching the observable
        variable_element = None
        for variable in module_component.findall('variable'):
            if variable.attrib['name'] == observable:
                variable_element = variable
                break      
        if variable_element is not None:
            variable_attributes = variable_element.attrib
            if 'public_interface' not in variable_attributes.keys():
                variable_attributes['public_interface']='out'
                copy_variable_attributes = copy.deepcopy(variable_attributes)
                copy_variable_attributes['public_interface']='in'
                ET.SubElement(model_component, 'variable', copy_variable_attributes)
            elif variable_attributes['public_interface'] == 'out':
                    copy_variable_attributes = copy.deepcopy(variable_attributes)
                    copy_variable_attributes['public_interface'] = 'in'
                    ET.SubElement(model_component, 'variable', copy_variable_attributes)
            else: # find the variable in the params component
                variable_attributes = {}
                for variable in params_component.findall('variable'):
                    if variable.attrib['name'] == observable:
                        variable_element = variable
                        break
                if variable_element is not None:
                    variable_attributes = variable_element.attrib                  
                    if variable_attributes['public_interface'] == 'out':
                        copy_variable_attributes = copy.deepcopy(variable_attributes)
                        copy_variable_attributes['public_interface'] = 'in'
                        # remove the initial value from the copy
                        if 'initial_value' in copy_variable_attributes.keys():
                            del copy_variable_attributes['initial_value']
                        ET.SubElement(model_component, 'variable', copy_variable_attributes)
                    else:
                        raise ValueError(f"Observable variable '{observable}' has an invalid public interface.")
                else:
                    raise ValueError(f"Observable variable '{observable}' not found")          

    module_element = import_comp(model_ET, module_name, module_name , module_file,module_ET)
    param_element= import_comp(model_ET, params_name, params_name , params_file,params_ET)    
    # map the variables in the module and main model
    map_components(model_ET, module_element, model_component)
    # map the variables in the parameters and main model
    map_components(model_ET, param_element, model_component)
    # map the variables in the module and parameters
    map_components(model_ET, module_element, param_element)
    
    # write to the CellML files
    module_filepath=Path(filepath).joinpath(module_file)
    params_filepath=Path(filepath).joinpath(params_file)
    model_filepath=Path(filepath).joinpath(main_file)
    write_cellmlV1(module_ET, module_filepath.resolve())
    write_cellmlV1(params_ET, params_filepath.resolve())
    write_cellmlV1(model_ET, model_filepath.resolve())

    return model_ET, module_ET, params_ET

def getRdfFile(filename):

  
    # Parse the XML
    tree = ET.parse(filename)
    root = tree.getroot()

    # Extract and print all namespaces (optional, for inspection)
    namespaces = dict([
        node for _, node in ET.iterparse(filename, events=['start-ns'])
    ])
    print("Detected namespaces:", namespaces)

    # register the namespaces
    for prefix, uri in namespaces.items():
        ET.register_namespace(prefix, uri)

    # Find <rdf:RDF> element(s)
    rdf_elements = root.findall(".//rdf:RDF", namespaces)
    output_filename = filename.replace('.cellml', '_rdf.xml')
    # Write to file
    with open(output_filename, 'w', encoding='utf-8') as out_file:
        out_file.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        for rdf in rdf_elements:
            out_file.write(ET.tostring(rdf, encoding='unicode'))

    print(f"RDF extracted and saved to {output_filename}")

def getAnother(causality):
    if causality=='effort':
        return 'flow'
    elif causality=='flow':
        return 'effort'
    else:
        raise ('causality is wrong')

def nxBG_CellML_elements(nxBGJson,full=False):
    CellML_elements={}
    CellML_elements['conservation_laws']=[]
    if not full:
        G = load_nxBG_json(nxBGJson)
        for node in G.nodes:
            if G.nodes[node]['a']=='BondElement':
                node_ID=node
                CellML_elements[node_ID]={}
                if G.nodes[node]['subClass']=='Se':
                    if 'modelParameter' in G.nodes[node]:
                        for param in G.nodes[node]['modelParameter'].keys():
                            if 'quantity' in G.nodes[node]['modelParameter'][param]['description']:
                                CellML_elements[node_ID]['quantity']=G.nodes[node]['modelParameter'][param]['propertyName']
                if 'modelState' in G.nodes[node]:
                    for state in G.nodes[node]['modelState']:
                        CellML_elements[node_ID]['quantity']=G.nodes[node]['modelState'][state]['propertyName']
                powerPorts=getPowerPorts(G,node)
                for powerPort in powerPorts:
                    CellML_elements[node_ID][powerPort]={}
                    causality=G.nodes[powerPort]['causality']
                    CellML_elements[node_ID][powerPort][causality]=G.nodes[powerPort][causality]['propertyName']
                    if 'expression' in G.nodes[powerPort][getAnother(causality)]:
                        eq={}
                        eq['yvar']=G.nodes[powerPort][getAnother(causality)]['propertyName']
                        eq['infix']=G.nodes[powerPort][getAnother(causality)]['expression']
                        CellML_elements['conservation_laws'].append(eq)
    return CellML_elements                       

if __name__ == "__main__": 
    
    nx_BG_file = './data/nx_BG_refined.json'
    G = load_nxBG_json(nx_BG_file)
    implementation_dict={'language':'CellML','version':1.1, 'filepath':'./data/',
                    'module_name':'SLC5_BG','module_file':'SLC5_BG.cellml',
                    'params_name':'SLC5_BG_param','params_file':'SLC5_BG_param.cellml',
                    'model_name':'SLC5_BG_run','model_file':'SLC5_BG_run.cellml',
                    'units_file':'SLC5_BG_units.cellml', 'voi':{'propertyName': 't', 'units': 'second'},
                    'observed_vars': {'v_r1','T'}}
    
    #model_ET, module_ET, params_ET=nxBG2CellMLV1(G, implementation_dict)

   
    getRdfFile('MacKenzie_1996.cellml')