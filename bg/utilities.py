import json
import numpy as np
import csv
import libsbml
from pathlib import Path
import warnings
#Get the directory where THIS script file lives
SCRIPT_DIR = Path(__file__).resolve().parent



def load_json(json_file: str, folder_path: str = '.././examples') -> dict:
    """
    Load the json file to a dictionary

    Parameters
    ----------
    json_file : str
        The file path of the json file
    folder_path : str
        The relative path of the folder where the json file is located
        By default, the json file will be loaded from the examples folder

    Returns
    -------
    comp_dict : dict
        The dictionary of the json file  

    """
    # check if the folder exists, if not, raise an error
    target_folder = Path(SCRIPT_DIR / folder_path).resolve()
    if not target_folder.exists():
        raise FileNotFoundError(f"The folder '{target_folder}' does not exist.")

    full_path = Path(target_folder / json_file)
    if not full_path.exists():
        raise FileNotFoundError(f"The file '{full_path}' does not exist.")

    with open(full_path, 'r') as f:
        comp_dict = json.load(f)

    return comp_dict

def save_json(data: dict, json_file: str, folder_path: str = '.././examples'):
    """
    Save the dictionary to a json file

    Parameters
    ----------
    data : dict
        The dictionary of the bond graph model
    json_file : str
        The name of the json file
    folder_path : str
        The relative path of the folder where the json file will be saved
        By default, the json file will be saved in the examples folder

    Returns
    -------
    None

    side effect
    ------------
    Save the dictionary to a json file

    """
    # create the folder if it does not exist
    target_folder = Path(SCRIPT_DIR / folder_path).resolve()
    if not target_folder.exists():
        target_folder.mkdir(parents=True, exist_ok=True)

    # check if the json file already exists, if yes, raise a warning
    full_path = Path(target_folder/ json_file)
    if full_path.exists():
        # Take user input to decide whether to overwrite the file or not
        user_input = input(f"The file '{full_path}' already exists. Do you want to overwrite it? (y/n): ")
        if user_input.lower() != 'y':
            print("Operation cancelled.")
            return
    with open(full_path, 'w') as f:
        json.dump(data, f,indent=4)

def load_matrix(matrix):
    """
    Load stoichiometric matrix from a csv file
    The ID is metamodel of the BG components.
    The components in the row take the potential as the input and the flow as the output, and we call them F components.
    The components in the column take the flow as the input and the potential as the output, and we call them E components.
    The csv file should have the following format (* means blank):
    *       *     *      fName fName
    *       *     *       fID   fID
    *       *     *      fPort fPort
    eName  eID  ePort      0    1 
    eName  eID  ePort      1    0
    
    Parameters
    ----------
    matrix : str
        The file path of the stoichiometric matrix
    
    Returns
    -------
    eName : list
        A list of E component (e_out, f_in) names
    eID : list
        A list of E component (e_out, f_in) IDs
    ePort : list
        A list of E component (e_out, f_in) port number
    fName : list
        A list of F component (e_in, f_out) names
    fID : list
        A list of F component (e_in, f_out) IDs
    fPort : list
        A list of F component (e_in, f_out) port number
    N : numpy.ndarray
        The stoichiometric matrix
    """
    startC=3
    N = []
    eName=[]
    eID=[]
    ePort=[]
    with open(matrix,'r') as f:
        reader = csv.reader(f,delimiter=',')
        line_count = 0
        for row in reader:
            if line_count ==0:
                fName=row[startC:]
                line_count += 1
            elif line_count ==1:
                fID=row[startC:]
                line_count += 1
            elif line_count == 2:
                fPort=row[startC:]
                line_count += 1
            else:
                N.append(row[startC:])
                eName.append(row[0])
                eID.append(row[1])
                ePort.append(row[2])
        f.close()
    
    return eName, eID, ePort, fName, fID, fPort, np.array(N).astype(int)

def load_matrix_domain(matrix):
    """
    Load stoichiometric matrix from a csv file
    The ID is metamodel of the BG components.
    The components in the row take the potential as the input and the flow as the output, and we call them F components.
    The components in the column take the flow as the input and the potential as the output, and we call them E components.
    The csv file should have the following format (* means blank):
    *      *       *     *      domain domain
    *      *       *     *      fName fName
    *      *       *     *       fID   fID
    *      *       *     *      fPort fPort
    domain  eName  eID  ePort      0    1 
    domain  eName  eID  ePort      1    0
    
    Parameters
    ----------
    matrix : str
        The file path of the stoichiometric matrix
    
    Returns
    -------
    eDomain : list
        A list of E component (e_out, f_in) domains
    eName : list
        A list of E component (e_out, f_in) names
    eID : list
        A list of E component (e_out, f_in) type IDs
    ePort : list
        A list of E component (e_out, f_in) port number
    fDomain : list
        A list of F component (e_in, f_out) domains
    fName : list
        A list of F component (e_in, f_out) names
    fID : list
        A list of F component (e_in, f_out) IDs
    fPort : list
        A list of F component (e_in, f_out) port number
    N : numpy.ndarray
        The stoichiometric matrix
    """
    startC=4
    N = []
    eDomain=[]
    eName=[]
    eID=[]
    ePort=[]
    with open(matrix,'r') as f:
        reader = csv.reader(f,delimiter=',')
        line_count = 0
        for row in reader:
            if line_count ==0:
                fDomain=row[startC:]
                line_count += 1
            elif line_count ==1:
                fName=row[startC:]
                line_count += 1
            elif line_count == 2:
                fID=row[startC:]
                line_count += 1
            elif line_count == 3:
                fPort=row[startC:]
                line_count += 1
            else:
                N.append(row[startC:])
                eDomain.append(row[0])
                eName.append(row[1])
                eID.append(row[2])
                ePort.append(row[3])
        f.close()
    
    return eDomain, eName, eID, ePort, fDomain, fName, fID, fPort, np.array(N).astype(int)

def infix_to_mathml(ode_var,infix, voi='',version='1.1'):
    """
    Convert the infix string to mathML string defined in CellML specification

    Parameters
    ----------
    ode_var : str
        The derivative variable name
    infix : str
        The infix string to be converted to mathML string
    voi : str, optional
        The variable of integration. The default is ''.
    version : str, optional
        The version of the CellML specification. The default is '1.1'.

    Returns
    -------
    str
        The mathML string defined in CellML specification

    """

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
    if '<cn type="integer">' in mathstr or '<cn type="real">' in mathstr or '<cn>' in mathstr:
        mathstr = mathstr.replace ('<cn type="integer">', '<cn cellml:units="dimensionless">')
        mathstr = mathstr.replace ('<cn type="real">', '<cn cellml:units="dimensionless">')
        mathstr = mathstr.replace ('<cn>', '<cn cellml:units="dimensionless">')
        # add left side of the equation       
    mathstr = preforumla + mathstr + postformula
    # add the cellml namespace to the mathML string
    if '<cn cellml:units="dimensionless">' in mathstr:
        mathstr = f'<apply xmlns:cellml="http://www.cellml.org/cellml/{version}#">' + mathstr + ' </apply>'
    else:
        mathstr = '<apply>' + mathstr + ' </apply>'
    return mathstr




def read_ParamCellML(sbml_file):
    """
    Read the parameter values from a CellML SBML file

    Parameters
    ----------
    sbml_file : str
        The file path of the CellML SBML file

    Returns
    -------
    param_dict : dict
        A dictionary of parameter names and their values

    """
    reader = libsbml.SBMLReader()
    document = reader.readSBML(sbml_file)
    model = document.getModel()
    param_dict = {}
    for parameter in model.getListOfParameters():
        param_dict[parameter.getId()] = parameter.getValue()
    return param_dict
    
        
    