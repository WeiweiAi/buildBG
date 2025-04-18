import copy
import numpy as np
import csv
from pathlib import Path
from utilities import *
from sympy import *
import pandas as pd
from scipy import integrate

defUnit=["ampere","becquerel","candela","celsius","coulomb","dimensionless","farad","gram","gray","henry",
    "hertz","joule","katal","kelvin","kilogram","liter","litre","lumen","lux","meter","metre","mole",
    "newton","ohm","pascal","radian","second","siemens","sievert","steradian","tesla","volt","watt","weber"]
params_common=['R','T','F']

def load_matrix(matrix):
    """
    Load stoichiometric matrix from a csv file
    The components are defined in BG_components.
    The ID of the components is the key of the dictionary BG_components.
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

def add_BG_components(bg_dict,bg_components_json, matrix, suffix=True):
    """
    Add the components from the stoichiometric matrix to the bond graph model
    The components are defined in BG_components.
    The ID of the components is the key of the dictionary BG_components.
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
    bg_dict : dict
        The dictionary of the bond graph model
        It can be empty or have some components
    bg_components_json : str
        The file path of the json file of the BG components definition
    matrix : str
        The file path of the stoichiometric matrix
    suffix : bool, optional
        If True, the symbols of the parameters, variables and state variables are updated with the component name as the suffix
        The default is True.
    
    Returns
    -------
    None

    side effect
    ------------
    Update the bg_dict dictionary with the components from the stoichiometric matrix
    """
    eName, eID, ePort, fName, fID, fPort,N=load_matrix(matrix)
    bg_components=load_json(bg_components_json)
    cNames=eName+fName
    cIDs=eID+fID
    # update the symbols of the parameters, variables and state variables for the components  
    for i in range(len(cNames)):
        cID=cIDs[i]
        if cID not in bg_components.keys():
            raise ValueError('The component type is not found in the BG_components')
        cIndex=cNames[i]# cIndex is the key for the bg_dict dictionary, make it unique; 
        if cIndex not in bg_dict.keys(): # only add the component if it is not in the dictionary
            bg_dict[cIndex]=copy.deepcopy(bg_components[cID])
            # Rename the parameters, variables and state variables with the component name as the suffix if suffix=True
            if suffix:
                if 'params' in bg_dict[cIndex].keys():
                    for param in bg_dict[cIndex]['params'].keys():
                        if param not in params_common:
                            bg_dict[cIndex]['params'][param]['symbol']=bg_dict[cIndex]['params'][param]['symbol']+ '_' + cIndex
                if 'vars' in bg_dict[cIndex].keys():
                    for var in bg_dict[cIndex]['vars'].keys():
                        bg_dict[cIndex]['vars'][var]['symbol']=bg_dict[cIndex]['vars'][var]['symbol']+ '_' + cIndex
        else:
            print('The component {} is already in the dictionary'.format(cIndex))

def add_BG_connections(bg_dict, matrix, direction):
    """
    Add the connections from the stoichiometric matrix to the bond graph model
    The components are defined in BG_components.
    The ID of the components is the key of the dictionary BG_components.
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
    bg_dict : dict
        The dictionary of the bond graph model
    matrix : str
        The file path of the stoichiometric matrix
    direction : str
        The direction of the flow
        'e2f': from E (0 node) to F (1 node)
        'f2e': from F (1 node) to E (0 node)
    
    Returns
    -------
    None

    side effect
    ------------
    Update the bg_dict dictionary with the connections from the stoichiometric matrix
    """
    eName, eID, ePort, fName, fID, fPort,N=load_matrix(matrix)
    # update the connections of the components
    for j in range(len(fName)):
        fIndex=fName[j]
        portF=fPort[j]
        if fIndex not in bg_dict.keys():
            raise ValueError('The component {} is not found in the BG dictionary'.format(fIndex))
        if portF not in bg_dict[fIndex]['ports'].keys():
            raise ValueError('The port {} is not found in the component {}'.format(portF,fIndex))
        for i in range(len(eName)):
            eIndex=eName[i]
            portE = ePort[i]
            if eIndex not in bg_dict.keys():
                raise ValueError('The component {} is not found in the BG dictionary'.format(eIndex))
            if portE not in bg_dict[eIndex]['ports'].keys():
                raise ValueError('The port {} is not found in the component {}'.format(portE,eIndex))
            if N[i,j]!=0:
                if bg_dict[fIndex]['ports'][portF]['type']=='Power' and bg_dict[eIndex]['ports'][portE]['type']=='Power':
                    bg_dict[fIndex]['ports'][portF]['causality']='e_in'
                    bg_dict[eIndex]['ports'][portE]['causality']='e_out'                   
                if direction == 'e2f':
                    if bg_dict[fIndex]['ports'][portF]['type']=='Signal':
                        bg_dict[fIndex]['ports'][portF]['causality']='s_in'
                    if bg_dict[eIndex]['ports'][portE]['type']=='Signal':
                        bg_dict[eIndex]['ports'][portE]['causality']='s_out'
                    if N[i,j]>0:
                        bg_dict[fIndex]['ports'][portF]['in']+=[[eIndex, portE, str(N[i,j])]]
                        bg_dict[eIndex]['ports'][portE]['out']+=[[fIndex, portF, str(N[i,j])]]
                    else:
                        bg_dict[fIndex]['ports'][portF]['out']+=[[eIndex, portE, str(-N[i,j])]]
                        bg_dict[eIndex]['ports'][portE]['in']+=[[fIndex, portF, str(-N[i,j])]]
                       
                elif direction == 'f2e':
                    if bg_dict[fIndex]['ports'][portF]['type']=='Signal':
                        bg_dict[fIndex]['ports'][portF]['causality']='s_out'
                    if bg_dict[eIndex]['ports'][portE]['type']=='Signal':
                        bg_dict[eIndex]['ports'][portE]['causality']='s_in'                                                                   
                    if N[i,j]>0:
                        bg_dict[fIndex]['ports'][portF]['out']+=[[eIndex, portE, str(N[i,j])]]
                        bg_dict[eIndex]['ports'][portE]['in']+=[[fIndex, portF, str(N[i,j])]]
                    else:
                        bg_dict[fIndex]['ports'][portF]['in']+=[[eIndex, portE, str(-N[i,j])]]
                        bg_dict[eIndex]['ports'][portE]['out']+=[[fIndex, portF, str(-N[i,j])]]
                else:
                    raise ValueError('The direction is not correct') 

def update_BG_eqn(bg_dict,voi):
    """
    Get the inputs for the components and update the constitutive relations of the components
    The inputs are the sum of the inputs from the ports
    The constitutive relations are updated with the real inputs

    Parameters
    ----------
    bg_dict : dict
        The dictionary of the bond graph model
    voi : dict
        The dictionary of the variables of integration
        The default is {
                "description": "Time",
                "units": "second",
                "symbol": "t"
        }.

    Returns
    -------
    None

    side effect
    ------------
    Update the constitutive relations of the components

    """
    def sum_bonds(sub_comps, causality):
        multiports=0
        for comp in sub_comps:
            cIndex=comp[0]
            portN=comp[1]
            stochoimetry=comp[2]       
            try:
                num_stochoimetry=int(stochoimetry)
            except:
                try:
                    num_stochoimetry=float(stochoimetry)
                except:
                    raise ValueError('The stochoimetry is not an integer or a float')
            if bg_dict[cIndex]['ports'][portN]['causality']=='e_out' and causality=='e_in':
                multiports+=num_stochoimetry*symbols(bg_dict[cIndex]['vars']['e_'+portN]['symbol'])   
            elif bg_dict[cIndex]['ports'][portN]['causality']=='e_in' and causality=='e_out':
                multiports+=num_stochoimetry*symbols(bg_dict[cIndex]['vars']['f_'+portN]['symbol'])
            elif bg_dict[cIndex]['ports'][portN]['causality']=='s_out': 
                pass
            elif bg_dict[cIndex]['ports'][portN]['causality']=='s_in':
                pass
            else:
                raise ValueError('The port number or variable causality is not correct')
                        
        return multiports
    
    def signal_in(sub_comps):
        signal=0
        for comp in sub_comps:
            cIndex=comp[0]
            portN=comp[1]
            stochoimetry=comp[2]       
            try:
                num_stochoimetry=int(stochoimetry)
            except:
                try:
                    num_stochoimetry=float(stochoimetry)
                except:
                    raise ValueError('The stochoimetry is not an integer or a float')
            if bg_dict[cIndex]['ports'][portN]['causality']=='e_out':
                signal+=num_stochoimetry*symbols(bg_dict[cIndex]['vars']['e_'+portN]['symbol'])   
            elif bg_dict[cIndex]['ports'][portN]['causality']=='e_in':
                signal+=num_stochoimetry*symbols(bg_dict[cIndex]['vars']['f_'+portN]['symbol'])
            elif bg_dict[cIndex]['ports'][portN]['causality']=='s_out':
                signal+=num_stochoimetry*symbols(bg_dict[cIndex]['vars']['s_'+portN]['symbol'])
            else:
                raise ValueError('The port number or variable causality is not correct')
                        
        return signal
    
    for key, comp in bg_dict.items():
        comp['conservation_eqs']=[]
        # get the inputs for the ports
        for port in comp['ports']:
            if comp['ports'][port]['type']=='Power':
                sub_comps_in=comp['ports'][port]['in']
                causality=comp['ports'][port]['causality']
                bonds_inward=sum_bonds(sub_comps_in,causality)         
                sub_comps_out=comp['ports'][port]['out']
                bonds_outward=sum_bonds(sub_comps_out,causality)
                if comp['ports'][port]['direction']=='in':
                    inputs_sum=bonds_inward-bonds_outward
                else:
                    inputs_sum=bonds_outward-bonds_inward
                 # replace the input variables with the real inputs
                if causality=='e_in':
                    comp['vars']['e_'+port]['expression']=ccode(inputs_sum)
                    comp['conservation_eqs']+= [[ comp['vars']['e_'+port]['symbol'],ccode(inputs_sum), '']]
                elif causality=='e_out':
                    comp['vars']['f_'+port]['expression']=ccode(inputs_sum)
                    comp['conservation_eqs']+= [[ comp['vars']['f_'+port]['symbol'],ccode(inputs_sum), '']]

            elif comp['ports'][port]['type']=='Signal':
                sub_comps_in=comp['ports'][port]['in']
                signal_in_=signal_in(sub_comps_in)
                comp['vars']['s_'+port]['expression']=ccode(signal_in_)
                comp['conservation_eqs']+= [[ comp['vars']['s_'+port]['symbol'],ccode(signal_in_), '']]
            else:
                raise ValueError('The port description is not correct')            
        comp['constitutive_eqs']=[]
        if 'vars' in comp.keys():
            var_dict_= copy.deepcopy(comp['vars'])
        if 'params' in comp.keys():
            var_dict_.update( copy.deepcopy(comp['params']))
        # update the constitutive relations template
        for str_expr in comp['constitutive_relations']:
            if 'ode' in str_expr:
                yvar=str_expr.split('ode(')[1].split(',')[0]
                #voi=str_expr.split('ode(')[1].split(',')[1].split(')')[0]
                other_terms=str_expr.split('ode(')[1].split(')')[1]
                expr=-sympify(other_terms) # assume the ode is on the left hand side
                for ikey,ivar in var_dict_.items():
                    if 'expression' in var_dict_[ikey].keys():
                        expr=expr.subs(symbols(ikey),sympify(var_dict_[ikey]['expression']))
                    else:
                        expr=expr.subs(symbols(ikey),symbols(var_dict_[ikey]['symbol']))
                comp['constitutive_eqs']+= [[comp['vars'][yvar]['symbol'], ccode(expr),  voi['symbol']]]
                if 'voi' not in comp['vars'].keys():
                    comp['vars']['voi']=voi
            else:
                expr=sympify(str_expr)
                for port in comp['ports']:
                    if comp['ports'][port]['causality']=='e_in' and symbols('f_'+port) in expr.free_symbols:
                       var_=solve(expr, symbols('f_'+port))
                       if len(var_)==0:
                            raise ValueError('The equation is not solvable')
                       else:
                            solved_expr=var_[0]
                       for ikey,ivar in var_dict_.items():
                            if 'expression' in var_dict_[ikey].keys():
                                solved_expr=solved_expr.subs(symbols(ikey),sympify(var_dict_[ikey]['expression']))
                            else:
                                solved_expr=solved_expr.subs(symbols(ikey),symbols(var_dict_[ikey]['symbol']))
                       comp['constitutive_eqs']+= [[comp['vars']['f_'+port]['symbol'], ccode(solved_expr), '']]
                       
                    if comp['ports'][port]['causality']=='e_out' and symbols('e_'+port) in expr.free_symbols:
                        var_=solve(expr, symbols('e_'+port))
                        if len(var_)==0:
                            raise ValueError('The equation is not solvable')
                        else:
                            solved_expr=var_[0]
                        for ikey,ivar in var_dict_.items():
                            if 'expression' in var_dict_[ikey].keys():
                                solved_expr=solved_expr.subs(symbols(ikey),symbols(var_dict_[ikey]['expression']))
                            else:
                                solved_expr=solved_expr.subs(symbols(ikey),symbols(var_dict_[ikey]['symbol']))
                        comp['constitutive_eqs']+= [[comp['vars']['e_'+port]['symbol'], ccode(solved_expr), '']]
                    
                    if comp['ports'][port]['causality']=='s_out' and symbols('s_'+port) in expr.free_symbols:
                        var_=solve(expr, symbols('s_'+port))
                        solved_expr=var_[0]
                        for ikey,ivar in var_dict_.items():
                            if 'expression' in var_dict_[ikey].keys():
                                solved_expr=solved_expr.subs(symbols(ikey),symbols(var_dict_[ikey]['expression']))
                            else:
                                solved_expr=solved_expr.subs(symbols(ikey),symbols(var_dict_[ikey]['symbol']))
                        comp['constitutive_eqs']+= [[comp['vars']['s_'+port]['symbol'], ccode(solved_expr), '']]                                  

def kinetic2BGparams(N_f,N_r,kf,kr,Kc,Nc,Ws):
    """
    Convert kinetic parameters to BG parameters for biochemical reactions
    The method is based on the thesis:
    Pan, Michael. A bond graph approach to integrative biophysical modelling.
    Diss. University of Melbourne, Parkville, Victoria, Australia, 2019.

    Parameters
    ----------
    N_f : numpy.ndarray
        The forward stoichiometry matrix
    N_r : numpy.ndarray
        The reverse stoichiometry matrix
    kf : 1d list
        The forward rate constants, the same order as the reactions in N_f  
    kr : 1d list
        The reverse rate constants, the same order as the reactions in N_r
    Kc : 1d list, the constraints, can be empty
    Nc : 2d list, can be empty, 
        the length of the Nc is the number of Kc,
        the length of the Nc[0] is the number of the species
    Ws : 1d list, the size is the number of species
           
    Returns
    -------
    kappa : 1d numpy.ndarray
        The reaction rate constants
    K : 1d numpy.ndarray
        The thermodynamic constants
    K_eq : numpy.ndarray
        The equilibrium constants
    diff_ : float
        The difference between the estimated and the input kinetic parameters
    zero_est : numpy.ndarray
        The estimated zero values of the detailed balance constraints

    """ 
    k_f=np.array([kf]).transpose()
    k_r=np.array([kr]).transpose()
    K_c=np.array([Kc]).transpose()
    N_c=np.array(Nc).transpose()
    W_s=np.array([Ws]).transpose()
    
    N_fT=np.transpose(N_f)
    N_rT=np.transpose(N_r)
    N = N_r - N_f
    num_cols = N_f.shape[1] # number of reactions, the same as the number of columns in N_f
    num_rows = N_f.shape[0] # number of species, the same as the number of rows in N_f
    I=np.identity(num_cols)
    N_cT=np.transpose(N_c)
    num_contraints = K_c.shape[0]
    zerofill=np.zeros((num_contraints,num_cols))
    K_eq = np.divide(k_f,k_r)
    if len(K_c)!=0:
        M=np.block([
            [I, N_fT],
            [I, N_rT],
            [zerofill, N_cT]
        ])
        k= np.block([
            [k_f],
            [k_r],
            [K_c]
        ]) 
        N_b =np.hstack([-N, N_c])
        K_contraints = np.block([
            [K_eq],
            [K_c]
        ])
    else:
        M=np.block([
            [I, N_fT],
            [I, N_rT]
        ])
        k= np.block([
            [k_f],
            [k_r]
        ])
        N_b = -N
        K_contraints = K_eq
    # construct W matrix 
    # the first nr elements are 1 for reactions, the last ns elements are the volume of species
    W=np.vstack([np.ones((num_cols,1)),W_s]) 

    # convert kinetic parameters to BG parameters
    lambdaW= np.exp(np.matmul(np.linalg.pinv(M),np.log(k)))
    lambda_ = np.divide(lambdaW,W)
    kappa=lambda_[:num_cols]
    K = lambda_[num_cols:]
    
    # check if the solution is valid
    N_rref, _ = Matrix(N).rref()
    zero_est = None
    R_mat = np.array(nsimplify(Matrix(-N), rational=True).nullspace())
    if R_mat.size>0:
        R_mat = np.transpose(np.array(R_mat).astype(np.float64))[0]
        zero_est = np.matmul(R_mat.T,K_eq)
    # Check that there is a detailed balance constraint
    if N_c.size>0:
        Z = np.array(nsimplify(Matrix(N_c), rational=True).nullspace()) #rational_nullspace(M, 2)
        if Z.size>0:
            Z = np.transpose(np.array(Z).astype(np.float64))[0]
            zero_est = np.matmul(Z.T,np.log(K_c))

    k_est = np.exp(np.matmul(M,np.log(lambdaW)))
    diff_ = np.sum(np.abs(np.divide(k_est - k,k)))

    return kappa[:,0], K[:,0], K_eq[:,0], diff_, zero_est,k_est

def update_BioBG_params(bg_dict, kappa, fName, K, eName):
    """
    Update the BG parameters of biochemitry components in the bg_dict
    Assume the parameters are kappa and K for the reactions and species components, respectively
    

    Parameters
    ----------
    bg_dict : dict
        The dictionary of the bond graph model
    kappa : numpy.ndarray
        The reaction rate constants
    fName : list
        A list of reaction component Re names
        The order of the names should be the same as the order of the kappa
    K : numpy.ndarray
        The thermodynamic constants
    eName : list
        A list of species component Ce names
        The order of the names should be the same as the order of the kappa

    Returns
    -------
    None

    side effect
    ------------
    Update the parameters of the biochemitry components
    """
    for i in range(len(fName)):
        if abs(kappa[i])<1e-2 or abs(kappa[i])>1e2:
            kappa_="{:.3e}".format(kappa[i])
        else:
            kappa_="{:.3f}".format(kappa[i])
        bg_dict[fName[i]]['params']['kappa']['value']=kappa_
    for i in range(len(eName)):
        if abs(K[i])<1e-2 or abs(K[i])>1e2:
            K_="{:.3e}".format(K[i])
        else:
            K_="{:.3f}".format(K[i])
        bg_dict[eName[i]]['params']['K']['value']=K_    

def update_BG_params(bg_dict, params):
    """
    Update the BG parameters in the bg_dict

    Parameters
    ----------
    bg_dict : dict
        The dictionary of the bond graph model
    params : a list of truple
        The parameters to update [(comp_name, param_name, param_val)] or [(param_name, param_val)]
    
    Returns
    -------
    None

    side effect
    ------------
    Update the parameters of the component

    """

    for i in range(len(params)):
        param=params[i]
        comp_name=None
        if len(param)==2:
            param_name=param[0]
            param_val=param[1]
        elif len(param)==3:
            param_name=param[1]
            param_val=param[2]
            comp_name=param[0]
        else:
            raise ValueError('The parameter format is not correct')
        if comp_name is None:
            for key, comp in bg_dict.items():
                if 'params' in comp.keys():
                    if param_name in comp['params'].keys():
                        bg_dict[key]['params'][param_name]['value']=param_val
        else:
            if comp_name not in bg_dict.keys():
                raise ValueError('The component {} is not found in the BG dictionary'.format(comp_name))
            if 'params' in bg_dict[comp_name].keys():
                bg_dict[comp_name]['params'][param_name]['value']=param_val

def write_params_csv( param_name,param_val, param_units, csv_file='params_BG.csv'):
    """
    Update the BG parameters in the bg_dict with the parameters

    Parameters
    ----------
    param_val : 1d numpy.ndarray
        The parameter value
    param_name : str array
        The parameter name
    param_units : str array
        The parameter units
    csv : str, optional
        The file path of the csv file to save the parameters
        The default is 'params_BG.csv'.

    Returns
    -------
    None

    side effect
    ------------
    Save the parameters to a csv file

    """
    
    with open(csv_file, mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(['Parameter', 'Value', 'Units'])
        for i in range(len(param_val)):
            writer.writerow([param_name[i], param_val[i], param_units[i]])
        file.close()

def build_AdjacencyMatrix_(bg_dict):
    """
    Build the adjacency matrix for the bond graph model
    Derive the sympy expression of the power on each bond in the model
    Derive the sympy expression of the power flowing into each component in the model

    Parameters
    ----------
    bg_dict : dict
        The dictionary of the bond graph model
        Requires the bond information for each component:
        comp['ports'][port]['in'] is a list of the inward connections
        comp['ports'][port]['out'] is a list of the outward connections
        comp['ports'][port]['causality'] define the direction of the effort variable of the port (e_in, e_out)
        comp['ports'][port]['direction'] is the positive direction of flow (in or out) relative to the component
        comp['vars']['e_'+port]['expression'] is the sympy expression of the potential on the port, optional
        comp['vars']['f_'+port]['expression'] is the sympy expression of the flow on the port, optional
        comp['vars']['e_'+port]['symbol'] is the sympy symbol of the potential on the port
        comp['vars']['f_'+port]['symbol'] is the sympy symbol of the flow on the port

    Returns
    -------
    A : numpy.ndarray
        The adjacency matrix, A[i][j]=1 if there is a connection from i to j
    comp_port:  list of str
        The list of the nodes (i.e.,'component,port#') 
    P_bond_expr : 2d list of sympy expression
        The sympy expression of power on each bond in the model
        P_bond_expr[i][j] is the sympy expression of power on the bond from i to j 
    P_comp_expr : dict
        The energy expression of each component in the model
        key is the component name
        P_comp_expr[key] is the sympy expression of power flowing into component key  
    """ 
    comp_port=[]
    for key, comp in bg_dict.items():
        if 'ports' in comp.keys():
            for port in comp['ports']:
                comp_port+=[key+','+port]

    A=np.zeros((len(comp_port),len(comp_port)))
    # create a dictionary for the power expression of each component
    P_comp_expr={}
    e_comp_port_expr=[] # save the list of effort expression of the ports
    f_comp_port_expr=[] # save the list of flow expression of the ports
    # create 2d list for the power expression
    P_bond_expr=[['' for i in range(len(comp_port))] for j in range(len(comp_port))]
    for key, comp in bg_dict.items():
        if 'ports' in comp.keys():
            P_sum=0
            for port in comp['ports']:
                e_str=comp['vars']['e_'+port]['expression'] if 'expression' in comp['vars']['e_'+port].keys() else comp['vars']['e_'+port]['symbol']
                f_str=comp['vars']['f_'+port]['expression'] if 'expression' in comp['vars']['f_'+port].keys() else comp['vars']['f_'+port]['symbol']
                P_symp=sympify(e_str)*sympify(f_str)
                e_comp_port_expr+=[sympify(e_str)]
                f_comp_port_expr+=[sympify(f_str)]
                if comp['ports'][port]['direction']=='in':
                    P_sum+=P_symp
                elif comp['ports'][port]['direction']=='out':
                    P_sum+=-P_symp
                else:
                    raise ValueError('The port direction is not correct')
                indexi=comp_port.index(key+','+port) # index of the port in the list B
                for i in range(len(comp['ports'][port]['in'])): # inwards: indexj to indexi
                    cIndex=comp['ports'][port]['in'][i][0]
                    portN=comp['ports'][port]['in'][i][1]
                    stoich=comp['ports'][port]['in'][i][2]
                    try:
                        num_stochoimetry=int(stoich)
                    except:
                        try:
                            num_stochoimetry=float(stoich)
                        except:
                            raise ValueError('The stochoimetry is not an integer or a float')                    
                    indexj=comp_port.index(cIndex+','+portN)
                    A[indexj][indexi]=1
                    if comp['ports'][port]['causality']=='e_in':
                        e_bond_str=num_stochoimetry*sympify(bg_dict[cIndex]['vars']['e_'+portN]['expression'] if 'expression' in bg_dict[cIndex]['vars']['e_'+portN].keys() else bg_dict[cIndex]['vars']['e_'+portN]['symbol'])
                        f_bond_str=sympify(f_str)
                    elif comp['ports'][port]['causality']=='e_out':
                        e_bond_str=sympify(e_str)
                        f_bond_str=num_stochoimetry*sympify(bg_dict[cIndex]['vars']['f_'+portN]['expression'] if 'expression' in bg_dict[cIndex]['vars']['f_'+portN].keys() else bg_dict[cIndex]['vars']['f_'+portN]['symbol'])
                    P_bond_expr[indexj][indexi]=e_bond_str*f_bond_str
                for i in range(len(comp['ports'][port]['out'])): # outwards: indexi to indexj
                    cIndex=comp['ports'][port]['out'][i][0]
                    portN=comp['ports'][port]['out'][i][1]
                    stoich=comp['ports'][port]['out'][i][2]
                    try:
                        num_stochoimetry=int(stoich)
                    except:
                        try:
                            num_stochoimetry=float(stoich)
                        except:
                            raise ValueError('The stochoimetry is not an integer or a float')
                    indexj=comp_port.index(cIndex+','+portN)
                    A[indexi][indexj]=1
                    if comp['ports'][port]['causality']=='e_in':
                        e_bond_str=num_stochoimetry*sympify(bg_dict[cIndex]['vars']['e_'+portN]['expression'] if 'expression' in bg_dict[cIndex]['vars']['e_'+portN].keys() else bg_dict[cIndex]['vars']['e_'+portN]['symbol'])
                        f_bond_str=sympify(f_str)
                    elif comp['ports'][port]['causality']=='e_out':
                        e_bond_str=sympify(e_str)
                        f_bond_str=num_stochoimetry*sympify(bg_dict[cIndex]['vars']['f_'+portN]['expression'] if 'expression' in bg_dict[cIndex]['vars']['f_'+portN].keys() else bg_dict[cIndex]['vars']['f_'+portN]['symbol'])
                    P_bond_expr[indexi][indexj]=e_bond_str*f_bond_str
            P_comp_expr[key]=P_sum      
    return A, comp_port, P_bond_expr, P_comp_expr, e_comp_port_expr, f_comp_port_expr 

def calc_energy(bg_json,result_csv):
    """
    Calculate the energy and activity of each bond in the model

    Parameters
    ----------
    bg_json : str, the fullpath to a json file
        The json file of the bond graph model

    result_csv : str, the fullpath to a csv file
        The csv file with the following format:
        t, var1, var2, ...
        Note that the first column should be the time variable
        Note that the variables should be the same as the symbols in the bond graph model


    Returns
    -------
    None

    side effect
    ------------
    save the power vector of each bond to a csv file 'bond_power.csv'
    save the power vector of each component to a csv file 'comp_power.csv'
    save the activity of each component to a csv file 'comp_activity.csv'    
    save the effort of each port to a csv file '_comp_port_e.csv'
    save the flow of each port to a csv file '_comp_port_f.csv'
    save the adjacency matrix to a csv file 'adjacency.csv'
    Save the energy and activity of each component to a json file 'activity.json'
    """ 
    bg_dict=load_json(bg_json)
    simResults_df= pd.read_csv(result_csv)
    # get the path of the csv file
    csv_path=Path(result_csv).parent
    csv_file_name=Path(result_csv).stem 
    csv_file_power_bond=csv_path/(csv_file_name+'_bond_power.csv')
    csv_file_power_comp=csv_path/(csv_file_name+'_comp_power.csv')
    csv_file_activity_comp=csv_path/(csv_file_name+'_comp_activity.csv')
    csv_file_e_comp_port=csv_path/(csv_file_name+'_comp_port_e.csv')
    csv_file_f_comp_port=csv_path/(csv_file_name+'_comp_port_f.csv')     
    csv_file_A=csv_path/(csv_file_name+'_adjacency.csv')  
    json_file=csv_path/(csv_file_name+'_activity.json')

    A, comp_port, P_bond_expr, P_comp_expr, e_comp_port_expr, f_comp_port_expr = build_AdjacencyMatrix_(bg_dict)
    df_A=pd.DataFrame(A)    
    df_A.to_csv(csv_file_A, index=False)
 
    E_comp_val=np.zeros(len(P_comp_expr))
    A_comp_val=np.zeros(len(P_comp_expr))
    AI_comp_val=np.zeros(len(P_comp_expr))    
    
    # create dataframes for the power vector of each bond and component
    df_power_bond=pd.DataFrame()
    df_e_comp_port=pd.DataFrame()
    df_f_comp_port=pd.DataFrame()    
    df_power_comp=pd.DataFrame(columns=bg_dict.keys())
    df_activity_comp=pd.DataFrame(columns=bg_dict.keys())
    df_power_bond['t']=simResults_df['t']
    df_e_comp_port['t']=simResults_df['t']
    df_f_comp_port['t']=simResults_df['t']
    df_power_comp['t']=simResults_df['t']
    df_activity_comp['t']=simResults_df['t']


    bond_ij=[]# save the list of bond (i,j)
    E_bond=[] # save the list of  bond energy
    A_bond=[] # save the list of  bond activity
    P_bond_expr_=[]# save the list of bond power expression
    P_comp_expr_=[]# save the list of component power expression
    e_comp_port_expr_=[]# save the list of effort expression of the ports
    f_comp_port_expr_=[]# save the list of flow expression of the ports
    
    for i in range(len(comp_port)):
        # calculate the energy and activity of each port
        list_vars=list(e_comp_port_expr[i].free_symbols)
        list_vars_str=[str(var) for var in list_vars]
        e_symp_func=lambdify(list_vars,e_comp_port_expr[i],'numpy')
        df_e_comp_port[i]=e_symp_func(*[simResults_df[var] for var in list_vars_str])
        list_vars=list(f_comp_port_expr[i].free_symbols)
        list_vars_str=[str(var) for var in list_vars]
        f_symp_func=lambdify(list_vars,f_comp_port_expr[i],'numpy')
        df_f_comp_port[i]=f_symp_func(*[simResults_df[var] for var in list_vars_str])
        e_comp_port_expr_+=[comp_port[i]+': '+str(e_comp_port_expr[i])]
        f_comp_port_expr_+=[comp_port[i]+': '+str(f_comp_port_expr[i])]
        # calculate the energy and activity of each bond
        for j in range(len(comp_port)):
            if A[i][j]==1:
                # translate the sympy expression to python function based on free symbols
                list_vars=list(P_bond_expr[i][j].free_symbols)
                list_vars_str=[str(var) for var in list_vars]
                P_symp_func=lambdify(list_vars,P_bond_expr[i][j],'numpy')
                P_bond_vec=P_symp_func(*[simResults_df[var] for var in list_vars_str])
                E_bond+=[np.trapz(P_bond_vec,simResults_df['t'])]
                A_bond+=[np.trapz(np.abs(P_bond_vec),simResults_df['t'])]
                df_power_bond[f'{i}_{j}']=P_bond_vec 
                bond_ij+=[(i,j)] 
                P_bond_expr_+=[comp_port[i]+'-->'+comp_port[j]+':   '+str(P_bond_expr[i][j])]

    # calculate the energy and activity of each component              
    A_total=0
    i=0
    for key in P_comp_expr.keys():
        list_vars=list(P_comp_expr[key].free_symbols)
        list_vars_str=[str(var) for var in list_vars]
        P_symp_func=lambdify(list_vars,P_comp_expr[key],'numpy')
        P_comp_vec=P_symp_func(*[simResults_df[var] for var in list_vars_str])
        E_comp_val[i]=np.trapz(P_comp_vec,simResults_df['t'])
        A_comp_val[i]=np.trapz(np.abs(P_comp_vec),simResults_df['t'])
        A_total+=A_comp_val[i]
        df_power_comp[key]=P_comp_vec
        df_activity_comp[key]=integrate.cumulative_trapezoid(np.abs(P_comp_vec),x=simResults_df['t'],initial=0)
        P_comp_expr_+=[key+':   '+str(P_comp_expr[key])]
        i+=1
           
    AI_comp_val=A_comp_val/A_total*100
    
    dict_activity={'Componentl list':list(bg_dict.keys()),'Energy': E_comp_val.tolist(), 'Activity': A_comp_val.tolist(), 'Activity Index': AI_comp_val.tolist(),
                   'Bond Energy': E_bond, 'Bond Activity': A_bond, 'Bond_ij': bond_ij, 'Comp_port': comp_port,
                   'P_bond_expr': P_bond_expr_, 'P_comp_expr': P_comp_expr_, 'e_comp_port': e_comp_port_expr_, 'f_comp_port': f_comp_port_expr_}
    
    # save the results to csv and json files
    df_power_bond.to_csv(csv_file_power_bond, index=False)
    df_power_comp.to_csv(csv_file_power_comp, index=False)
    df_e_comp_port.to_csv(csv_file_e_comp_port, index=False)
    df_f_comp_port.to_csv(csv_file_f_comp_port, index=False)
    df_activity_comp.to_csv(csv_file_activity_comp, index=False)    
    save_json(dict_activity, json_file)   

    return 


if __name__ == "__main__": 
    
    # specify the csv files for the stoichiometric matrix
    file_path='./data/'
    fmatrix=file_path+'SLC5_f.csv'
    rmatrix=file_path+'SLC5_r.csv'
    bg_components_json='BG_components.json'        
    
    # build the bond graph model
    bg_dict={}       
    direction = 'e2f'
    add_BG_components(bg_dict,bg_components_json, fmatrix)
    add_BG_connections(bg_dict, fmatrix, direction)
    direction = 'f2e'
    add_BG_connections(bg_dict, rmatrix, direction)
    
    # update the equations of the bond graph model
    voi={'description': 'Time', 'units': 'second', 'symbol': 't'}
    update_BG_eqn(bg_dict,voi)
    
    fmatrix=file_path+'SLC5_f_chem.csv'
    rmatrix=file_path+'SLC5_r_chem.csv' 
    # update the BG parameters for the biochemical reactions
    eName, eID, ePort, fName, fID, fPort,N_f=load_matrix(fmatrix)
    eName, eID, ePort, fName, fID, fPort,N_r=load_matrix(rmatrix)
    AVO=6.022e23
    V1=1e-3
    k_12=8e4*V1*V1
    k_23=1e5*V1
    k_34=50
    k_45=800
    k_56=10
    k_61=5
    k_25=0.3
    k_21=500
    k_32=20
    k_43=50
    k_65=50*V1*V1
    k_16=35
    k_52=k_12*k_25*k_56*k_61/(k_21*k_65*k_16)
    k_54=k_23*k_34*k_45*k_52/(k_32*k_43*k_25)

    balance1=k_12*k_23*k_34*k_45*k_56*k_61/(k_21*k_32*k_43*k_54*k_65*k_16)
    balance2=k_12*k_25*k_56*k_61/(k_21*k_52*k_65*k_16)
    balance3=k_23*k_34*k_45*k_52/(k_32*k_43*k_54*k_25)
    print(balance1,balance2,balance3)

    kf=[k_12, k_23, k_34, k_45, k_56, k_61, k_25]
    kr=[k_21, k_32, k_43, k_54, k_65, k_16, k_52 ]
    K_c=[]
    N_c=[]

    V_E=1
    V_o=8.5e5
    V_i=8.5e5
    Ws=[V_i,V_o,V_i,V_o,V_E,V_E,V_E,V_E,V_E,V_E]
    kappa, K, K_eq, diff_,  zero_est,k_est= kinetic2BGparams(N_f,N_r,kf,kr,K_c,N_c,Ws)
    print(diff_,zero_est,K_eq,k_est)
    update_BioBG_params(bg_dict, kappa, fName, K, eName)
    update_BG_params(bg_dict,[('T',310)])
    # combine K and kappa to a single array
    param_val=np.concatenate((kappa,K))
    param_name=fName+eName
    param_units=['$fmol/s$']*len(fName)+['$fmol^{-1}$']*len(eName)
    write_params_csv(param_name,param_val,param_units,csv_file=file_path+'params_BG.csv')
    # save the bond graph model to a json file   
    json_file=file_path+'SLC5_BG.json'    
    save_json(bg_dict, json_file)   

    # calculate the energy and activity of the bond graph model
    file_path=r'C:\Users\wai484\temp\Energy-based-System-Analysis\models\\'
   # calc_energy(file_path+'SLC5_BG.json',file_path+'report_task_SGLT1_BGEA.csv')




        

        
