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


def kinetic2BGparams(N_f,N_r,kf,kr,K_c,N_c,Ws):
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
    kf : numpy.ndarray
        The forward rate constants,
        a column vector with the same number of rows as the number of reactions
        and the same order as the reactions in N_f  
    kr : numpy.ndarray
        The reverse rate constants,
        a column vector with the same number of rows as the number of reactions
        and the same order as the reactions in N_r
    K_c : numpy.ndarray
        The constraints vector,
        a column vector
    N_c : numpy.ndarray
        The constraints matrix,
        the columns of N_c is the same as the number of the K_c
        the rows of N_c is the same as the number of the species
    Ws : numpy.ndarray
        The volume vector, the size is the number of species ns       

    Returns
    -------
    kappa : numpy.ndarray
        The reaction rate constants
    K : numpy.ndarray
        The thermodynamic constants
    K_eq : numpy.ndarray
        The equilibrium constants
    diff_ : float
        The difference between the estimated and the input kinetic parameters
    zero_est : numpy.ndarray
        The estimated zero values of the detailed balance constraints

    """ 
    N_fT=np.transpose(N_f)
    N_rT=np.transpose(N_r)
    N = N_r - N_f
    num_cols = N_f.shape[1] # number of reactions, the same as the number of columns in N_f
    num_rows = N_f.shape[0] # number of species, the same as the number of rows in N_f
    I=np.identity(num_cols)
    N_cT=np.transpose(N_c)
    num_contraints = K_c.shape[0]
    zerofill=np.zeros((num_contraints,num_cols))
    K_eq = np.divide(kf,kr)
    if len(K_c)!=0:
        M=np.block([
            [I, N_fT],
            [I, N_rT],
            [zerofill, N_cT]
        ])
        k= np.block([
            [kf],
            [kr],
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
            [kf],
            [kr]
        ])
        N_b = -N
        K_contraints = K_eq
    # construct W matrix 
    # the first nr elements are 1 for reactions, the last ns elements are the volume of species
    W=np.vstack([np.ones((num_cols,1)),Ws]) 

    # convert kinetic parameters to BG parameters
    lambdaW= np.exp(np.matmul(np.linalg.pinv(M),np.log(k)))
    lambda_ = np.divide(lambdaW,W)
    kappa=lambda_[:num_cols]
    K = lambda_[num_cols:]
    
    # check if the solution is valid
    N_rref, _ = Matrix(N).rref()
    zero_est = None
    R_mat = np.array(nsimplify(Matrix(N), rational=True).nullspace())
    if R_mat.size>0:
        R_mat = np.transpose(np.array(R_mat).astype(np.float64))[0]
        zero_est = np.matmul(R_mat.T,K_eq)
    # Check that there is a detailed balance constraint
    Z = nsimplify(Matrix(N_b), rational=True).nullspace() #rational_nullspace(M, 2)
    if Z:
        Z = np.transpose(np.array(Z).astype(np.float64))[0]
        zero_est = np.matmul(Z.T,np.log(K_contraints))

    k_est = np.exp(np.matmul(M,np.log(lambdaW)))
    diff_ = np.sum(np.abs(np.divide(k_est - k,k)))

    return kappa, K, K_eq, diff_, zero_est

def update_BG_params(comp_dict, kappa, RName, K, CName, csv_file='params_BG.csv'):
    """
    Update the BG parameters in the comp_dict with the parameters

    Parameters
    ----------
    comp_dict : dict
        The dictionary of the bond graph model
    kappa : numpy.ndarra
        The reaction rate constants
    RName : list
        A list of reaction component Re names
        The order of the names should be the same as the order of the kappa
    K : numpy.ndarray
        The thermodynamic constants
    CName : list
        A list of species component Ce names
        The order of the names should be the same as the order of the kappa
    csv : str, optional
        The file path of the csv file to save the parameters
        The default is 'params_BG.csv'.

    Returns
    -------
    None

    side effect
    ------------
    Update the parameters of the components
    Save the parameters to a csv file

    """
    def toLatexStr (var_name):
        sub_str_1=var_name.split('_')[0]
        if len(var_name.split('_'))>1:
            sub_str_2=var_name.split('_')[1]
            latex_str='$\\'+sub_str_1+'_{'+sub_str_2+'}$'
        else:
            sub_str_2=''
            latex_str='$\\'+sub_str_1+'$'

        return latex_str
    
    with open(csv_file, mode='w') as file:
        writer = csv.writer(file)
        writer.writerow(['Parameter', 'Value'])
        for i in range(len(RName)):
            if abs(kappa[i][0])<1e-2 or abs(kappa[i][0])>1e2:
                kappa_="{:.3e}".format(kappa[i][0])
            else:
                kappa_="{:.3f}".format(kappa[i][0])
            comp_dict[RName[i]]['params']['kappa']['value']=kappa_
            var_name=comp_dict[RName[i]]['params']['kappa']['symbol']
            latex_str=toLatexStr(var_name)
            writer.writerow([latex_str, kappa_])
        for i in range(len(CName)):
            if abs(K[i][0])<1e-2 or abs(K[i][0])>1e2:
                K_="{:.3e}".format(K[i][0])
            else:
                K_="{:.3f}".format(K[i][0])
            comp_dict[CName[i]]['params']['K']['value']=K_
            var_name=comp_dict[CName[i]]['params']['K']['symbol']
            latex_str=toLatexStr(var_name)
            writer.writerow([latex_str, K_])

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
    
    CName, CType, CPort, RName, RType, RPort,N_f=load_matrix(file_path+fmatrix)
    CName, CType, CPort, RName, RType, RPort,N_r=load_matrix(file_path+rmatrix)

    V=1
    V_o=90
    h=0.726;g=12.1;c=1113;d=90.3;a=500000*V_o;b=a*9.5;f=3000*V_o;e=12.8459*f
    kf=np.array([[h, c, a, e]]).transpose()
    kr=np.array([[g, d, b, f]]).transpose()
    K_c=np.array([[1]]).transpose()
    N_c=np.array([[1,-1,0,0,0,0]]).transpose()
    K_c=np.array([[]]).transpose()
    N_c=np.array([[]]).transpose()
    V_i=0.09
    V_o=0.09
    V_E=1
    Ws=np.array([[V_i,V_o,V_E,V_E,V_E,V_E]]).transpose()

    kappa, K, K_eq, diff_,  zero_est= kinetic2BGparams(N_f,N_r,kf,kr,K_c,N_c,Ws)
    update_BG_params(comp_dict, kappa, RName, K, CName, csv_file='params_BG.csv')

    json_file=Path(fmatrix).with_suffix('.json')    
    with open(json_file, 'w') as f:
        json.dump(comp_dict, f,indent=4)   

        

        
