import pandas as pd
import copy
import numpy as np
import csv
from pathlib import Path
import json
from sympy import *

defUnit=["ampere","becquerel","candela","celsius","coulomb","dimensionless","farad","gram","gray","henry",
    "hertz","joule","katal","kelvin","kilogram","liter","litre","lumen","lux","meter","metre","mole",
    "newton","ohm","pascal","radian","second","siemens","sievert","steradian","tesla","volt","watt","weber"]
params_common=['R','T','F']

def load_matrix(matrix):
    """
    Load stoichiometric matrices from csv files
    The components are defined in BG_components.
    The type of the components is the key of the dictionary BG_components.
    The components in the row take the potential as the input and the flow as the output, and we call them R components.
    The components in the column take the flow as the input and the potential as the output, and we call them C components.
    The csv file should have the following format (* means blank):
    *       *     *      RName RName
    *       *     *      RType RType
    *       *     *      RPort RPort
    CName CType CPort      0    1 
    CName CType CPort      1    0
    
    Parameters
    ----------
    matrix : str
        The file path of the forward stoichiometric matrix
    
    Returns
    -------
    CName : list
        A list of C component (e_out, f_in) names
    CType : list
        A list of C component (e_out, f_in) types
    CPort : list
        A list of C component (e_out, f_in) port number
    RName : list
        A list of R component (e_in, f_out) names
    RType : list
        A list of R component (e_in, f_out) types
    RPort : list
        A list of R component (e_in, f_out) port number
    N : numpy.ndarray
        The stoichiometric matrix
    """
    startC=3
    N = []
    CName=[]
    CType=[]
    CPort=[]
    with open(matrix,'r') as f:
        reader = csv.reader(f,delimiter=',')
        line_count = 0
        for row in reader:
            if line_count ==0:
                RName=row[startC:]
                line_count += 1
            elif line_count ==1:
                RType=row[startC:]
                line_count += 1
            elif line_count == 2:
                RPort=row[startC:]
                line_count += 1
            else:
                N.append(row[startC:])
                CName.append(row[0])
                CType.append(row[1])
                CPort.append(row[2])
        f.close()
    
    return CName, CType, CPort, RName, RType, RPort, np.array(N).astype(int)

def build_BG_Dict(matrix,bg_components,comp_dict, direction,file_path='./'):
    """
    Build the dictionary for the components and connections from the stoichiometric matrix
    Parameters
    ----------
    matrix : str
        The file path of the stoichiometric matrix
        The csv file should have the following format (* means blank):
        *       *     *      RName RName
        *       *     *      RType RType
        *       *     *      RPort RPort
        CName CType CPort      0    1
        CName CType CPort      1    0
    bg_components : dict
        The dictionary of the bond graph components (predefined)
    comp_dict : dict
        The dictionary of the bond graph model
        Can be empty or have some components already
    direction : str
        The direction of the flow
        'C2R': from C (0 node) to R (1 node)
        'R2C': from R (1 node) to C (0 node)
    file_path : str, optional
        The file path of the csv files
        The default is './'.
    Returns
    -------
    None

    side effect
    ------------
    Update the comp_dict dictionary
    """
    CName, CType, CPort, RName, RType, RPort,N=load_matrix(file_path+matrix)
    CNames=CName+RName
    CTypes=CType+RType
    # update the symbols of the parameters, variables and state variables for the components  
    for i in range(len(CNames)):
        CType=CTypes[i]
        CIndex=CNames[i]# CName is the key for the dictionary, make it unique; 
        if CType in bg_components.keys() and CIndex not in comp_dict.keys():
            comp_dict[CIndex]=copy.deepcopy(bg_components[CType])
            comp_dict[CIndex]['type']=CType
            comp_dict[CIndex]['constitutive_relations_sym']=[]
            # Rename the parameters, variables and state variables with the component name as the suffix
            for param in comp_dict[CIndex]['params'].keys():
                if param not in params_common:
                    comp_dict[CIndex]['params'][param]['symbol']=comp_dict[CIndex]['params'][param]['symbol']+ '_' + CIndex
            for var in comp_dict[CIndex]['vars'].keys():
                comp_dict[CIndex]['vars'][var]['symbol']=comp_dict[CIndex]['vars'][var]['symbol']+ '_' + CIndex
            if 'state_vars' in bg_components[CType].keys():                  
                for state_var in bg_components[CType]['state_vars'].keys():
                    comp_dict[CIndex]['state_vars'][state_var]['symbol']=comp_dict[CIndex]['state_vars'][state_var]['symbol']+ '_' + CIndex                 
        elif CType not in bg_components.keys():
            raise ValueError('The component type is not found in the BG_components')
        else:
            pass
    # update the connections of the components
    for j in range(len(RName)):
        RIndex=RName[j]
        portR=RPort[j]
        if portR not in ['0','1']:
            raise ValueError('The port number is not correct')
        for i in range(len(CName)):
            CIndex=CName[i]
            portC = CPort[i]
            if portC not in ['0','1']:
                raise ValueError('The port number is not correct')
            if N[i,j]!=0: 
                comp_dict[RIndex]['vars']['e_'+portR]['IO']='in'
                comp_dict[RIndex]['vars']['f_'+portR]['IO']='out'
                comp_dict[RIndex]['ports'][portR]['type']='e_in'
                comp_dict[CIndex]['vars']['e_'+portC]['IO']='out'
                comp_dict[CIndex]['vars']['f_'+portC]['IO']='in'
                comp_dict[CIndex]['ports'][portC]['type']='e_out'
                if direction == 'C2R':
                    if N[i,j]>0:
                        comp_dict[RIndex]['ports'][portR]['in']+=[[CIndex, portC, str(N[i,j])]]
                        comp_dict[CIndex]['ports'][portC]['out']+=[[RIndex, portR, str(N[i,j])]]
                    else:
                        comp_dict[RIndex]['ports'][portR]['out']+=[[CIndex, portC, str(N[i,j])]]
                        comp_dict[CIndex]['ports'][portC]['in']+=[[RIndex, portR, str(N[i,j])]]  
                elif direction == 'R2C':                                                                   
                    if N[i,j]>0:
                        comp_dict[RIndex]['ports'][portR]['out']+=[[CIndex, portC, str(N[i,j])]]
                        comp_dict[CIndex]['ports'][portC]['in']+=[[RIndex, portR, str(N[i,j])]]
                    else:
                        comp_dict[RIndex]['ports'][portR]['in']+=[[CIndex, portC, str(N[i,j])]]
                        comp_dict[CIndex]['ports'][portC]['out']+=[[RIndex, portR, str(N[i,j])]]
                else:
                    raise ValueError('The direction is not correct') 

def update_BG_eqn(comp_dict):
    """
    Get the inputs for the components and update the constitutive relations of the components
    The inputs are the sum of the inputs from the ports
    The constitutive relations are updated with the real inputs

    Parameters
    ----------
    comp_dict : dict
        The dictionary of the bond graph model

    Returns
    -------
    None

    side effect
    ------------
    Update the constitutive relations of the components

    """
    def get_inputs(sub_comps):
        inputs=0
        for comp in sub_comps:
            CIndex=comp[0]
            portN=comp[1]
            stochoimetry=comp[2]            
            try:
                num_stochoimetry=int(stochoimetry)
            except:
                try:
                    num_stochoimetry=float(stochoimetry)
                except:
                    raise ValueError('The stochoimetry is not an integer or a float')
            if comp_dict[CIndex]['ports'][portN]['type']=='e_out':
                inputs+=num_stochoimetry*symbols(comp_dict[CIndex]['vars']['e_'+portN]['symbol'])    
            elif comp_dict[CIndex]['ports'][portN]['type']=='e_in':
                inputs+=num_stochoimetry*symbols(comp_dict[CIndex]['vars']['f_'+portN]['symbol'])
            else:
                raise ValueError('The port number or variable type is not correct')
                        
        return inputs
    
    for key, comp in comp_dict.items():
         
        if 'vars' in comp.keys():
            var_dict_= copy.deepcopy(comp['vars'])
        if 'state_vars' in comp.keys():
            var_dict_.update( copy.deepcopy(comp['state_vars']))
        if 'params' in comp.keys():
            var_dict_.update( copy.deepcopy(comp['params']))
        # get the inputs for the ports
        for port in comp['ports']:
            sub_comps_in=comp['ports'][port]['in']
            inputs_inward=get_inputs(sub_comps_in)         
            sub_comps_out=comp['ports'][port]['out']
            inputs_outward=get_inputs(sub_comps_out)
            if comp['ports'][port]['direction']=='in':
                inputs_sum=inputs_inward-inputs_outward
            else:
                inputs_sum=inputs_outward-inputs_inward
             # replace the input variables with the real inputs
            if comp['ports'][port]['type']=='e_in':
                comp['vars']['e_'+port]['expression']=ccode(inputs_sum)
            elif comp['ports'][port]['type']=='e_out':
                comp['vars']['f_'+port]['expression']=ccode(inputs_sum)
        # update the constitutive relations template
        for str_expr in comp['constitutive_relations']:
            if 'ode' in str_expr:
                yvar=str_expr.split('ode(')[1].split(',')[0]
                voi=str_expr.split('ode(')[1].split(',')[1].split(')')[0]
                other_terms=str_expr.split('ode(')[1].split(')')[1]
                expr=-sympify(other_terms) # assume the ode is on the left hand side
                # replace the input variables with the real inputs
                for port in comp['ports']:
                    if comp['ports'][port]['type']=='e_in' and 'expression' in comp['vars']['e_'+port].keys():
                        expr=expr.subs(symbols('e_'+port),sympify(comp['vars']['e_'+port]['expression']))
                    elif comp['ports'][port]['type']=='e_out' and 'expression' in comp['vars']['f_'+port].keys():
                        expr=expr.subs(symbols('f_'+port),sympify(comp['vars']['f_'+port]['expression']))
                
                for ikey,ivar in var_dict_.items():
                    expr=expr.subs(symbols(ikey),symbols(var_dict_[ikey]['symbol']))

                comp['constitutive_relations_sym']+= [[ccode(expr), comp['state_vars'][yvar]['symbol'], comp['voi'][voi]['symbol']]]
            else:
                expr=sympify(str_expr)
                for key_var,var in comp['vars'].items():
                    if 'IO' in comp['vars'][key_var] and comp['vars'][key_var]['IO']=='out' and symbols(key_var) in expr.free_symbols:
                        var_=solve(expr, symbols(key_var))
                        solved_expr=var_[0]
                        for port in comp['ports']:
                            if (comp['ports'][port]['type']=='e_in' ) and 'expression' in comp['vars']['e_'+port].keys():                                
                                solved_expr=solved_expr.subs(symbols('e_'+port),sympify(comp['vars']['e_'+port]['expression']))
                            if (comp['ports'][port]['type']=='e_out' ) and 'expression' in comp['vars']['f_'+port].keys():
                                solved_expr=solved_expr.subs(symbols('f_'+port),sympify(comp['vars']['f_'+port]['expression']))
                               
                        for ikey,ivar in var_dict_.items():
                            solved_expr=solved_expr.subs(symbols(ikey),symbols(var_dict_[ikey]['symbol']))
                        comp['constitutive_relations_sym']+= [[ccode(solved_expr), comp['vars'][key_var]['symbol'], '']]

if __name__ == "__main__": 

    file_path='./data/'
    fmatrix='SLC2_f_1.csv'
    rmatrix='SLC2_r_1.csv'
    bg_components_json='BG_components.json'
    with open(bg_components_json) as f:
        bg_components = json.load(f)
    
    comp_dict={}
    direction = 'C2R'
    build_BG_Dict(fmatrix,bg_components,comp_dict,direction,file_path=file_path)
    direction = 'R2C'
    build_BG_Dict(rmatrix,bg_components,comp_dict,direction,file_path=file_path)
    
    update_BG_eqn(comp_dict)
    # dump the comp_dict to a json file, which has the same name and path as the csv file (using Path)
    json_file=Path(fmatrix).with_suffix('.json')    
    with open(json_file, 'w') as f:
        json.dump(comp_dict, f,indent=4)   

        

        
