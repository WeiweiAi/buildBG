import numpy as np
from sympy import *
import libsbml 
from utilities import *
import xml.etree.ElementTree as ET
from xml.dom import minidom
import copy

defUnit=["ampere","becquerel","candela","celsius","coulomb","dimensionless","farad","gram","gray","henry",
    "hertz","joule","katal","kelvin","kilogram","liter","litre","lumen","lux","meter","metre","mole",
    "newton","ohm","pascal","radian","second","siemens","sievert","steradian","tesla","volt","watt","weber"]
params_common=['R','T','F']

CellMLV1_namespaces = {
        'cellml': "http://www.cellml.org/cellml/1.1#",  # CellML namespace
        'xlink': "http://www.w3.org/1999/xlink",  # XLink namespace
        'math': "http://www.w3.org/1998/Math/MathML"  # MathML namespace
    }

def infix_to_mathml(ode_var,infix, voi=''):
    if voi!='':
        preforumla = '<apply> <eq/> <apply> <diff/> <bvar> <ci>'+ voi + '</ci> </bvar> <ci>' + ode_var + '</ci> </apply> '
    else:
        preforumla = '<apply> <eq/> <ci>'+ ode_var + '</ci>'    
    postformula = ' </apply> '
    # replace log to ln in infix string
    infix = infix.replace('log', 'ln')
    # replace fabs to abs in infix string
    infix = infix.replace('fabs', 'abs')
    p = libsbml.parseL3Formula (infix)
    mathstr = libsbml.writeMathMLToString (p)
     # remove the <math> tags in the mathML string, and the namespace declaration will be added later according to the CellML specification
    mathstr = mathstr.replace ('<math xmlns="http://www.w3.org/1998/Math/MathML">', '')
    mathstr = mathstr.replace ('</math>', '')
    mathstr = mathstr.replace ('<?xml version="1.0" encoding="UTF-8"?>', '')
    # temporary solution to add cellml units for constant in the mathML string, replace <cn type="integer"> to <cn cellml:units="dimensionless">
    # check if <cn type="integer"> is in the mathML string
    if '<cn type="integer">' in mathstr or '<cn type="real">' in mathstr:
        mathstr = mathstr.replace ('<cn type="integer">', '<cn cellml:units="dimensionless">')
        mathstr = mathstr.replace ('<cn type="real">', '<cn cellml:units="dimensionless">')
        # add left side of the equation       
    mathstr = preforumla + mathstr + postformula
    # add the cellml namespace to the mathML string
    if '<cn cellml:units="dimensionless">' in mathstr:
        mathstr = '<apply xmlns:cellml="http://www.cellml.org/cellml/1.1#">' + mathstr + ' </apply>'
    else:
        mathstr = '<apply>' + mathstr + ' </apply>'
    return mathstr

def variable_attributes_IO(variable_attributes, IO_string):
    if IO_string == 'pub_in':
        variable_attributes['public_interface'] = 'in'
    elif IO_string == 'pub_out':
        variable_attributes['public_interface'] = 'out'
    elif IO_string == 'priv_in':
        variable_attributes['private_interface'] = 'in'
    elif IO_string == 'priv_out':
        variable_attributes['private_interface'] = 'out'
    elif IO_string == '':
        pass
    

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
    
    return ET.Element('model', model_attrs)

def CellML_param_ET(bg_dict,model_ET):
    """Add parameters to a CellML V1.x model as ET. ElementTree

    Parameters
    ----------
    bg_dict : dict
        The dictionary of the bond graph model
    model_ET : xml.etree.ElementTree
        The model to which the parameters will be added

    Returns
    -------
    param_units_Set : set
        A set of units used in the parameters

    Side effects
    ------------
    The model_ET is modified in place

    """
    # register the CellML namespace    
    ET.register_namespace('cellml', CellMLV1_namespaces['cellml'])
    model_name=model_ET.attrib['name']   
    param_units_Set=set()     
    # Add component and param variables
    component = ET.SubElement(model_ET, 'component', {'name': model_name})
    params_set=set()
    for comp in bg_dict:
        for param in bg_dict[comp]['params']:
            param_units_Set.add(bg_dict[comp]['params'][param]['units'])
            if bg_dict[comp]['params'][param]['symbol'] not in params_set: # components may share parameters such as R, T, F
                variable_attributes = {
                    'name': bg_dict[comp]['params'][param]['symbol'],
                    'units': bg_dict[comp]['params'][param]['units']
                }
                if 'value' in bg_dict[comp]['params'][param]:
                    variable_attributes['initial_value'] = str(bg_dict[comp]['params'][param]['value'])
                if 'IO' in bg_dict[comp]['params'][param]:
                    IO_string = bg_dict[comp]['params'][param]['IO']
                    variable_attributes_IO(variable_attributes, IO_string)
                else:
                    variable_attributes['public_interface']='out' # default is public output                  
                ET.SubElement(component, 'variable', variable_attributes)
                params_set.add(bg_dict[comp]['params'][param]['symbol'])
            else:
                pass
    return param_units_Set   

def CellML_model_run(bg_dict,run_ET,observables):
    """Add observed variables to a CellML V1.x model as ET. ElementTree

    Parameters
    ----------
    bg_dict : dict
        The dictionary of the bond graph model
    run_ET : xml.etree.ElementTree
        The model to which the observed variables will be added

    observables : list of tuples of str
        The list of observed variables, each tuple contains
          the component name,
          the 'vars' key of the variable, 
          and 'IO': pub_in, pub_out, priv_in, priv_out, or none

    Returns
    -------
    var_units_Set : set
        A set of units used in the observed variables

    """
    # register the CellML namespace    
    ET.register_namespace('cellml', CellMLV1_namespaces['cellml'])
    model_name=run_ET.attrib['name']
    # Add component and variables
    component = ET.SubElement(run_ET, 'component', {'name': model_name})
    # Add variables and units
    units_set = set()
    for observable in observables:
        comp=observable[0]
        var=observable[1]
        IO_string=observable[2]
        if var in bg_dict[comp]['vars']:
            variable_attributes = {'name': bg_dict[comp]['vars'][var]['symbol'],
                    'units': bg_dict[comp]['vars'][var]['units']}
            variable_attributes_IO(variable_attributes, IO_string)
                
            units_set.add(bg_dict[comp]['vars'][var]['units'])
            ET.SubElement(component, 'variable', variable_attributes)
        if var in bg_dict[comp]['params']:
            variable_attributes = {'name': bg_dict[comp]['params'][var]['symbol'],
                    'units': bg_dict[comp]['params'][var]['units']}
            variable_attributes_IO(variable_attributes, IO_string)
            units_set.add(bg_dict[comp]['params'][var]['units'])
            ET.SubElement(component, 'variable', variable_attributes)
    return units_set    


def CellML_model_ET(bg_dict,model_ET):
    """Add model to a CellML V1.x model as ET. ElementTree

    Parameters
    ----------
    bg_dict : dict
        The dictionary of the bond graph model
    model_ET : xml.etree.ElementTree
        The model to which the model will be added

    Returns
    -------
    var_units_Set : set
        A set of units used in the model

    Side effects
    ------------
    The model_ET is modified in place

    """
    # register the CellML namespace    
    ET.register_namespace('cellml', CellMLV1_namespaces['cellml'])
    model_name=model_ET.attrib['name']
    # Add component and variables
    component = ET.SubElement(model_ET, 'component', {'name': model_name})
    # Create the MathML element with the correct namespace
    mathml_element = ET.Element('math', {'xmlns': CellMLV1_namespaces['math']})
    param_set = set()
    var_set = set()
    units_set = set()
    param_attrs = []
    var_attrs = []
    for comp in bg_dict:
        if 'params' not in bg_dict[comp]:
            bg_dict[comp]['params'] = {}
        for param in bg_dict[comp]['params']:
            if bg_dict[comp]['params'][param]['symbol'] not in param_set:
                variable_attributes = {
                    'name': bg_dict[comp]['params'][param]['symbol'],
                    'units': bg_dict[comp]['params'][param]['units']
                }
                if 'IO' in bg_dict[comp]['params'][param]:
                    IO_string = bg_dict[comp]['params'][param]['IO']
                    variable_attributes_IO(variable_attributes, IO_string)
                else:
                    variable_attributes['public_interface']='in' # default is public input
                units_set.add(bg_dict[comp]['params'][param]['units'])
                param_attrs.append(variable_attributes)
                param_set.add(bg_dict[comp]['params'][param]['symbol'])
            else:
                pass
        if 'vars' not in bg_dict[comp]:
            bg_dict[comp]['vars'] = {}
        for var in bg_dict[comp]['vars']:
            if 'expression' in bg_dict[comp]['vars'][var]:
                pass # skip variables that are defined by expressions
            else:
                if bg_dict[comp]['vars'][var]['symbol'] not in var_set:
                    var_set.add(bg_dict[comp]['vars'][var]['symbol'])
                    if 'initial value' in bg_dict[comp]['vars'][var]:
                        if isinstance(bg_dict[comp]['vars'][var]['initial value'], str):
                            initial_value_str = bg_dict[comp]['params'][bg_dict[comp]['vars'][var]['initial value']]['symbol']
                        else:
                            initial_value_str = str(bg_dict[comp]['vars'][var]['initial value'])
                        variable_attributes = {
                            'name': bg_dict[comp]['vars'][var]['symbol'],
                            'units': bg_dict[comp]['vars'][var]['units'],
                            'initial_value': initial_value_str
                        }
                    else:
                        variable_attributes = {
                            'name': bg_dict[comp]['vars'][var]['symbol'],
                            'units': bg_dict[comp]['vars'][var]['units']
                        }
                    units_set.add(bg_dict[comp]['vars'][var]['units'])
                    if 'IO' in bg_dict[comp]['vars'][var]:
                        IO_string = bg_dict[comp]['vars'][var]['IO']
                        variable_attributes_IO(variable_attributes, IO_string)                
                    var_attrs.append(variable_attributes)
        if 'constitutive_eqs' not in bg_dict[comp]:
            bg_dict[comp]['constitutive_eqs'] = []
        for equation in bg_dict[comp]['constitutive_eqs']:
            mmathml_string = infix_to_mathml(equation[0], equation[1], equation[2])
        # Parse the MathML string
            try:
                # Convert the generated math string into an XML element
                math_content = ET.fromstring(mmathml_string)
                if len(math_content) > 0:
                    mathml_element.append(math_content[0])  # Append the first child, which is the actual content
            except ET.ParseError as e:
                print(f"Error parsing MathML: {e}")
    # Add parameters
    for param_attr in param_attrs:
        ET.SubElement(component, 'variable', param_attr)
    # Add variables
    for var_attr in var_attrs:
        ET.SubElement(component, 'variable', var_attr)
    # Append the <math> element to the component
    component.append(mathml_element)
    return units_set

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
    units_import = ET.SubElement(model_ET, 'import')
    units_import.set('xlink:href', units_file)
    for units_name in units_Set:
        if units_name not in defUnit:
            units_import_i = ET.SubElement(units_import, 'units', {'name': units_name, 'units_ref': units_name})

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
            model_import = ET.SubElement(model_ET, 'import')
            model_import.set('xlink:href', model_file)
            import_component= ET.SubElement(model_import, 'component', {'name': comp_name_new, 'component_ref': comp_name})
            return component_element  
    if import_component is None:
        print('The component is not found in the imported model')
        return None

def map_components(model, component1, component2):
    """Map the variables of two components in a CellML V1.x model as ET. ElementTree

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
                    if  'public_interface' in variable_1.attrib:
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

def CellML_model_param_combine(model_ET,param_ET, run_ET, inforun):
    """
    Combine model_ET, param_ET and run_ET to create a CellML V1.x model as ET. ElementTree
    The mapping based on the name of the variables and also the units of the variables

    Parameters
    ----------
    model_ET : xml.etree.ElementTree
        The model containing the model variables and equations
    param_ET : xml.etree.ElementTree
        The parameters model containing the model parameters
    run_ET : xml.etree.ElementTree
        The run model containing the observed variables
    inforun : dict
        The dictionary of the run model information
        'component_name': component name of the model,
        'component_name_new': component name in the new model,
        'model_file': model file name including the path
        'params_name': component name of the parameters
        'params_name_new': component name in the new model
        'params_file': parameters file name including the path

    Returns
    -------
    None  

    """
    component_element = import_comp(run_ET, inforun['component_name_new'], inforun['component_name'], inforun['model_file'],model_ET)
    param_element= import_comp(run_ET, inforun['params_name_new'], inforun['params_name'], inforun['params_file'],param_ET)    
    run_element=run_ET.find('component') # assume that the run model has only one component  
    # map the variables in the model and run model
    map_components(run_ET, component_element, run_element)
    # map the variables in the parameters and run model
    map_components(run_ET, param_element, run_element)
    # map the variables in the model and parameters
    map_components(run_ET, component_element, param_element)

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

if __name__ == "__main__": 

     # save the bond graph model to a json file
    file_path='./data/'
    bg_dict=load_json(file_path+'SLC5_BG.json')
    implemention={'language':'CellML','version':1.1,
                    'model_name':'SLC5_BG','model_file':'SLC5_BG.cellml',
                    'params_model_name':'SLC5_BG_param','params_file':'SLC5_BG_param.cellml',
                    'model_sim':'SLC5_BG_run','sim_file':'SLC5_BG_run.cellml',
                    'observed_vars': [['r1','f0','pub_in']],}
    model_dict={'SLC5_BG':{'description': 'SLC5_BG model','implemention': implemention, 'bg_dict': bg_dict}} 
    
    observed_vars=[('r1','f_0','pub_in')]

    units_file='./units.cellml'
    param_ET=create_cellmlV1_rootET('SLC5_BG_param', cellml_prefix=True)
    model_ET=create_cellmlV1_rootET('SLC5_BG', cellml_prefix=False)
    run_ET=create_cellmlV1_rootET('SLC5_BG_run', cellml_prefix=True)
    param_units_Set=CellML_param_ET(bg_dict,param_ET)
    var_units_Set=CellML_model_ET(bg_dict,model_ET)
    run_units_Set=CellML_model_run(bg_dict,run_ET,observed_vars)
    units_import(model_ET, var_units_Set,units_file)
    units_import(param_ET, param_units_Set,units_file)
    units_import(run_ET, run_units_Set,units_file)
    inforun={'component_name': 'SLC5_BG', 'component_name_new': 'SLC5_BG', 'model_file': 'SLC5_BG.cellml',
             'params_name': 'SLC5_BG_param', 'params_name_new': 'SLC5_BG_param', 'params_file': 'SLC5_BG_param.cellml'}
    CellML_model_param_combine(model_ET,param_ET, run_ET, inforun)
    write_cellmlV1(param_ET,file_path+'SLC5_BG_param.cellml')
    write_cellmlV1(model_ET,file_path+'SLC5_BG.cellml')
    write_cellmlV1(run_ET,file_path+'SLC5_BG_run.cellml')
    write_cellmlV1(read_cellmlV1(file_path+'SLC5_BG.cellml'),file_path+'SLC5_BG_test.cellml')
    
    
  #  to_cellmlV1_params(comp_dict, model_name='params_BG',model_file='params_BG.txt')
  #  to_cellmlV1_models(comp_dict, model_name='GLUT2_BG',model_file='GLUT2_BG.txt',params_file='params_BG.cellml')