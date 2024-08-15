import numpy as np
from sympy import *
import libsbml 
from utilities import *
import xml.etree.ElementTree as ET
from xml.dom import minidom

""" Misc. Constants of CellML construction """
MATH_HEADER = '<math xmlns="http://www.w3.org/1998/Math/MathML" xmlns:cellml="http://www.cellml.org/cellml/2.0#">\n'
MATH_FOOTER = '</math>\n'

defUnit=["ampere","becquerel","candela","celsius","coulomb","dimensionless","farad","gram","gray","henry",
    "hertz","joule","katal","kelvin","kilogram","liter","litre","lumen","lux","meter","metre","mole",
    "newton","ohm","pascal","radian","second","siemens","sievert","steradian","tesla","volt","watt","weber"]
params_common=['R','T','F']

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
    mathstr = mathstr.replace ('<cn type="integer">', '<cn cellml:units="dimensionless">')
    # add left side of the equation       
    mathstr = preforumla + mathstr + postformula
    return mathstr

def json2CellMLV1_param(json_file, model_name, component_name):
    """Convert a JSON file to a CellML V1.x model as ET. ElementTree with only parameters

    Parameters
    ----------
    json_file : str
        The file path of the JSON file
    model_name : str
        The name of the model
    component_name : str
        The name of the component

    Returns
    -------
    model : xml.etree.ElementTree
        The model as an ElementTree    
    """

    comp_dict=load_json(json_file)
    # Define namespaces with the 'cellml' prefix explicitly included
    namespaces = {
        '': "http://www.cellml.org/cellml/1.1#",  # Default namespace
        'xlink': "http://www.w3.org/1999/xlink",
    }
    # Register namespaces
    for prefix, uri in namespaces.items():
        ET.register_namespace(prefix, uri)

    # Create the root model element with namespaces
    model = ET.Element(ET.QName(namespaces[''], 'model'), {
        'name': model_name
    })

    # Add units import
    units_import = ET.SubElement(model, ET.QName(namespaces[''], 'import'), {ET.QName(namespaces['xlink'], 'href'): './units.cellml'}) 
    param_units_Set=set()     
    # Add component and param variables
    component = ET.SubElement(model, ET.QName(namespaces[''], 'component'), {'name': component_name})
    params_set=set()
    for comp in comp_dict:
        for param in comp_dict[comp]['params']:
            param_units_Set.add(comp_dict[comp]['params'][param]['units'])
            if comp_dict[comp]['params'][param]['symbol'] not in params_set:
                variable_attributes = {
                    'name': comp_dict[comp]['params'][param]['symbol'],
                    'units': comp_dict[comp]['params'][param]['units']
                }
                if 'value' in comp_dict[comp]['params'][param]:
                    variable_attributes['initial_value'] = str(comp_dict[comp]['params'][param]['value'])
                if 'IO' in comp_dict[comp]['params'][param]:
                    variable_attributes['public_interface'] = comp_dict[comp]['params'][param]['IO']
                else:
                    variable_attributes['public_interface'] = 'out' # default is output
                ET.SubElement(component, ET.QName(namespaces[''], 'variable'), variable_attributes)
                params_set.add(comp_dict[comp]['params'][param]['symbol'])
            else:
                pass

    for units_name in param_units_Set:
        if units_name not in defUnit:
            units_import_i = ET.SubElement(units_import, ET.QName(namespaces[''], 'units'), {'name': units_name, 'units_ref': units_name})  

    return model

def json2CellMLV1_model(json_file, model_name, component_name):
    """Convert a JSON file to a CellML V1.x model as ET. ElementTree with parameters, 
       variables, state variables, and mathematical equations

    Parameters
    ----------
    json_file : str
        The file path of the JSON file
    model_name : str
        The name of the model
    component_name : str
        The name of the component

    Returns
    -------
    model : xml.etree.ElementTree
        The model as an ElementTree
    """
    comp_dict=load_json(json_file)
    # Define namespaces with the 'cellml' prefix explicitly included
    namespaces = {
        '': "http://www.cellml.org/cellml/1.1#",  # Default namespace
        'xlink': "http://www.w3.org/1999/xlink",
    }
    # Register namespaces
    for prefix, uri in namespaces.items():
        ET.register_namespace(prefix, uri)

    # Create the root model element with namespaces
    model = ET.Element(ET.QName(namespaces[''], 'model'), {
        'name': model_name
    })

    # Add units import
    units_import = ET.SubElement(model, ET.QName(namespaces[''], 'import'), {ET.QName(namespaces['xlink'], 'href'): './units.cellml'})
    # Collect all units used in the model and add them to the import
    units_set = set()    
    # Add component and variables
    component = ET.SubElement(model, ET.QName(namespaces[''], 'component'), {'name': component_name})

    # Create the <math> element
    math_ns = "http://www.w3.org/1998/Math/MathML"
    math_element = ET.Element(ET.QName(math_ns, 'math'))
    math_element.set('xmlns', math_ns)
    # Collect parameters, variables, and state variables 
    param_set = set()
    var_set = set()
    param_attrs = []
    var_attrs = []
    state_var_attrs = []
    for comp in comp_dict:
        if 'params' not in comp_dict[comp]:
            comp_dict[comp]['params'] = {}
        for param in comp_dict[comp]['params']:
            if comp_dict[comp]['params'][param]['symbol'] not in param_set:
                variable_attributes = {
                    'name': comp_dict[comp]['params'][param]['symbol'],
                    'units': comp_dict[comp]['params'][param]['units']
                }
                if 'IO' in comp_dict[comp]['params'][param]:
                    variable_attributes['public_interface'] = comp_dict[comp]['params'][param]['IO']
                else:
                    variable_attributes['public_interface'] = 'in' # default is input
                units_set.add(comp_dict[comp]['params'][param]['units'])
                param_attrs.append(variable_attributes)
                param_set.add(comp_dict[comp]['params'][param]['symbol'])
            else:
                pass
        if 'vars' not in comp_dict[comp]:
            comp_dict[comp]['vars'] = {}
        for var in comp_dict[comp]['vars']:
            if 'expression' in comp_dict[comp]['vars'][var]:
                pass # skip variables that are defined by expressions
            else:
                if comp_dict[comp]['vars'][var]['symbol'] not in var_set:
                    var_set.add(comp_dict[comp]['vars'][var]['symbol'])
                    variable_attributes = {
                        'name': comp_dict[comp]['vars'][var]['symbol'],
                        'units': comp_dict[comp]['vars'][var]['units']
                    }
                    units_set.add(comp_dict[comp]['vars'][var]['units'])
                    if 'IO' in comp_dict[comp]['vars'][var]:
                        variable_attributes['public_interface'] = comp_dict[comp]['vars'][var]['IO']                
                    var_attrs.append(variable_attributes)
                    
        if 'state_vars' in comp_dict[comp].keys():
            for state_var in comp_dict[comp]['state_vars']:
                q_init = comp_dict[comp]['state_vars'][state_var]['value']
                if isinstance(q_init, str):
                    q_init_str = comp_dict[comp]['params'][q_init]['symbol']
                else:
                    q_init_str = str(q_init)
                variable_attributes = {
                    'name': comp_dict[comp]['state_vars'][state_var]['symbol'],
                    'units': comp_dict[comp]['state_vars'][state_var]['units'],
                    'initial_value': q_init_str
                }
                units_set.add(comp_dict[comp]['state_vars'][state_var]['units'])
                if 'IO' in comp_dict[comp]['state_vars'][state_var]:
                    variable_attributes['public_interface'] = comp_dict[comp]['state_vars'][state_var]['IO']
                else:
                    pass # default is internal
                state_var_attrs.append(variable_attributes)
        if 'constitutive_eqs' not in comp_dict[comp]:
            comp_dict[comp]['constitutive_eqs'] = []
        for equation in comp_dict[comp]['constitutive_eqs']:
            mmathml_string = infix_to_mathml(equation[0], equation[1], equation[2])
        # Parse the MathML string
            mathml_element = ET.fromstring(mmathml_string)
            math_element.append(mathml_element)
        #for equation in comp_dict[comp]['conservation_eqs']:
        #    mmathml_string = infix_to_mathml(equation[0], equation[1], equation[2])
        # Parse the MathML string
        #    mathml_element = ET.fromstring(mmathml_string)
        #    math_element.append(mathml_element)
    # Add units import
    for units_name in units_set:
        if units_name not in defUnit:
            units_import_i = ET.SubElement(units_import, ET.QName(namespaces[''], 'units'), {'name': units_name, 'units_ref': units_name})

    # Add parameters
    for param_attr in param_attrs:
        ET.SubElement(component, ET.QName(namespaces[''], 'variable'), param_attr)
    # Add variables
    for var_attr in var_attrs:
        ET.SubElement(component, ET.QName(namespaces[''], 'variable'), var_attr)
    # Add state variables
    for state_var_attr in state_var_attrs:
        ET.SubElement(component, ET.QName(namespaces[''], 'variable'), state_var_attr)

     # Append the <math> element to the component
    component.append(math_element)

    return model



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
    rough_string = ET.tostring(model, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    pretty_string = reparsed.toprettyxml(indent="    ")

    # Write the pretty-printed XML to file
    with open(model_file, 'w') as output_file:
        output_file.write(pretty_string)    


if __name__ == "__main__": 

    param_model=json2CellMLV1_param('./data/SLC2_BG.json', 'test_model_param', 'param_component')
    test_model= json2CellMLV1_model('./data/SLC2_BG.json', 'test_model', 'component_test')
    write_cellmlV1(param_model,'test_model_param.cellml')
    write_cellmlV1(test_model,'test_model.cellml')
    
  #  to_cellmlV1_params(comp_dict, model_name='params_BG',model_file='params_BG.txt')
  #  to_cellmlV1_models(comp_dict, model_name='GLUT2_BG',model_file='GLUT2_BG.txt',params_file='params_BG.cellml')