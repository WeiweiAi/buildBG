import copy
import numpy as np
import csv
from pathlib import Path
from utilities import *
from sympy import *

defUnit=["ampere","becquerel","candela","celsius","coulomb","dimensionless","farad","gram","gray","henry",
    "hertz","joule","katal","kelvin","kilogram","liter","litre","lumen","lux","meter","metre","mole",
    "newton","ohm","pascal","radian","second","siemens","sievert","steradian","tesla","volt","watt","weber"]
params_common=['R','T','F']

def load_matrix(matrix):
    """
    Load stoichiometric matrices from csv files
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
        The file path of the forward stoichiometric matrix
    
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

def build_BG_Dict(matrix,bg_components,bg_dict, direction,file_path='./'):
    """
    Build the dictionary for the components and connections from the stoichiometric matrix
    Parameters
    ----------
    matrix : str
        The file path of the stoichiometric matrix
        The csv file should have the following format (* means blank):
         *       *     *      fName fName
         *       *     *       fID   fID
         *       *     *      fPort fPort
        eName   eID   ePort      0    1 
        eName   eID   ePort      1    0
    bg_components : dict
        The dictionary of the bond graph components (predefined)
    bg_dict : dict
        The dictionary of the bond graph model
        Can be empty or have some components already
    direction : str
        The direction of the flow
        'e2f': from E (0 node) to F (1 node)
        'f2e': from F (1 node) to E (0 node)
    file_path : str, optional
        The file path of the csv files
        The default is './'.
    Returns
    -------
    None

    side effect
    ------------
    Update the bg_dict dictionary
    """
    eName, eID, ePort, fName, fID, fPort,N=load_matrix(file_path+matrix)
    cNames=eName+fName
    cIDs=eID+fID
    # update the symbols of the parameters, variables and state variables for the components  
    for i in range(len(cNames)):
        cID=cIDs[i]
        cIndex=cNames[i]# CName is the key for the dictionary, make it unique; 
        if cID in bg_components.keys() and cIndex not in bg_dict.keys():
            bg_dict[cIndex]=copy.deepcopy(bg_components[cID])
            bg_dict[cIndex]['constitutive_eqs']=[]
            bg_dict[cIndex]['conservation_eqs']=[]
            # Rename the parameters, variables and state variables with the component name as the suffix
            for param in bg_dict[cIndex]['params'].keys():
                if param not in params_common:
                    bg_dict[cIndex]['params'][param]['symbol']=bg_dict[cIndex]['params'][param]['symbol']+ '_' + cIndex
            for var in bg_dict[cIndex]['vars'].keys():
                bg_dict[cIndex]['vars'][var]['symbol']=bg_dict[cIndex]['vars'][var]['symbol']+ '_' + cIndex
            if 'state_vars' in bg_components[cID].keys():                  
                for state_var in bg_components[cID]['state_vars'].keys():
                    bg_dict[cIndex]['state_vars'][state_var]['symbol']=bg_dict[cIndex]['state_vars'][state_var]['symbol']+ '_' + cIndex            
        elif cID not in bg_components.keys():
            raise ValueError('The component type is not found in the BG_components')
        else:
            pass
    # update the connections of the components
    for j in range(len(fName)):
        fIndex=fName[j]
        portF=fPort[j]
        if 's' in portF: # signal port number starts with s
            connectionF='signals'
        else:
            connectionF='ports'
        for i in range(len(eName)):
            eIndex=eName[i]
            portE = ePort[i]
            if 's' in portE:
                connectionE='signals'
            else:
                connectionE='ports'
            if N[i,j]!=0:
                if connectionF=='ports' and connectionE=='ports':
                    bg_dict[fIndex][connectionF][portF]['type']='e_in'
                    bg_dict[eIndex][connectionE][portE]['type']='e_out'
                if connectionF=='signals':
                   bg_dict[fIndex][connectionF][portF]['type']='s_in'
                if connectionE=='signals':
                   bg_dict[eIndex][connectionE][portE]['type']='s_out'  # signal from e to f                    
                if direction == 'e2f':
                    if N[i,j]>0:
                        bg_dict[fIndex][connectionF][portF]['in']+=[[eIndex, portE, str(N[i,j])]]
                        bg_dict[eIndex][connectionE][portE]['out']+=[[fIndex, portF, str(N[i,j])]]
                    else:
                        bg_dict[fIndex][connectionF][portF]['out']+=[[eIndex, portE, str(N[i,j])]]
                        bg_dict[eIndex][connectionE][portE]['in']+=[[fIndex, portF, str(N[i,j])]]
                       
                elif direction == 'f2e':                                                                   
                    if N[i,j]>0:
                        bg_dict[fIndex][connectionF][portF]['out']+=[[eIndex, portE, str(N[i,j])]]
                        bg_dict[eIndex][connectionE][portE]['in']+=[[fIndex, portF, str(N[i,j])]]
                    else:
                        bg_dict[fIndex][connectionF][portF]['in']+=[[eIndex, portE, str(N[i,j])]]
                        bg_dict[eIndex][connectionE][portE]['out']+=[[fIndex, portF, str(N[i,j])]]
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
    def sum_bonds(sub_comps):
        multiports=0
        for comp in sub_comps:
            cIndex=comp[0]
            portN=comp[1]
            stochoimetry=comp[2]
            if 's' in portN:
                connectionN='signals'
            else:
                connectionN='ports'          
            try:
                num_stochoimetry=int(stochoimetry)
            except:
                try:
                    num_stochoimetry=float(stochoimetry)
                except:
                    raise ValueError('The stochoimetry is not an integer or a float')
            if bg_dict[cIndex][connectionN][portN]['type']=='e_out':
                multiports+=num_stochoimetry*symbols(bg_dict[cIndex]['vars']['e_'+portN]['symbol'])   
            elif bg_dict[cIndex][connectionN][portN]['type']=='e_in':
                multiports+=num_stochoimetry*symbols(bg_dict[cIndex]['vars']['f_'+portN]['symbol'])
            elif bg_dict[cIndex][connectionN][portN]['type']=='s_out':
                pass
            elif bg_dict[cIndex][connectionN][portN]['type']=='s_in':
                pass
            else:
                raise ValueError('The port number or variable type is not correct')
                        
        return multiports
    
    def signal_in(sub_comps):
        signal=0
        for comp in sub_comps:
            cIndex=comp[0]
            portN=comp[1]
            stochoimetry=comp[2]
            if 's' in portN:
                connectionN='signals'
            else:
                connectionN='ports'          
            try:
                num_stochoimetry=int(stochoimetry)
            except:
                try:
                    num_stochoimetry=float(stochoimetry)
                except:
                    raise ValueError('The stochoimetry is not an integer or a float')
            if bg_dict[cIndex][connectionN][portN]['type']=='e_out':
                signal+=num_stochoimetry*symbols(bg_dict[cIndex]['vars']['e_'+portN]['symbol'])   
            elif bg_dict[cIndex][connectionN][portN]['type']=='e_in':
                signal+=num_stochoimetry*symbols(bg_dict[cIndex]['vars']['f_'+portN]['symbol'])
            elif bg_dict[cIndex][connectionN][portN]['type']=='s_out':
                signal+=num_stochoimetry*symbols(bg_dict[cIndex]['vars']['s_'+portN]['symbol'])
            else:
                raise ValueError('The port number or variable type is not correct')
                        
        return signal
    
    for key, comp in bg_dict.items():
        # get the inputs for the ports
        for port in comp['ports']:
            sub_comps_in=comp['ports'][port]['in']
            bonds_inward=sum_bonds(sub_comps_in)         
            sub_comps_out=comp['ports'][port]['out']
            bonds_outward=sum_bonds(sub_comps_out)
            if comp['ports'][port]['direction']=='in':
                inputs_sum=bonds_inward-bonds_outward
            else:
                inputs_sum=bonds_outward-bonds_inward
             # replace the input variables with the real inputs
            if comp['ports'][port]['type']=='e_in':
                comp['vars']['e_'+port]['expression']=ccode(inputs_sum)
                comp['conservation_eqs']+= [[ comp['vars']['e_'+port]['symbol'],ccode(inputs_sum), '']]
            elif comp['ports'][port]['type']=='e_out':
                comp['vars']['f_'+port]['expression']=ccode(inputs_sum)
                comp['conservation_eqs']+= [[ comp['vars']['f_'+port]['symbol'],ccode(inputs_sum), '']]
        
        if 'vars' in comp.keys():
            var_dict_= copy.deepcopy(comp['vars'])
        if 'state_vars' in comp.keys():
            var_dict_.update( copy.deepcopy(comp['state_vars']))
        if 'params' in comp.keys():
            var_dict_.update( copy.deepcopy(comp['params']))

        for signal in comp['signals']:
            sub_comps_in=comp['signals'][signal]['in']
            signal_in_=signal_in(sub_comps_in)
            comp['vars']['s_'+signal]['expression']=ccode(signal_in_)
            comp['constitutive_eqs']+= [[ comp['vars']['s_'+signal]['symbol'],ccode(signal_in_), '']]
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
                comp['constitutive_eqs']+= [[comp['state_vars'][yvar]['symbol'], ccode(expr),  voi['symbol']]]
            else:
                expr=sympify(str_expr)
                for port in comp['ports']:
                    if comp['ports'][port]['type']=='e_in' and symbols('f_'+port) in expr.free_symbols:
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
                       
                    if comp['ports'][port]['type']=='e_out' and symbols('e_'+port) in expr.free_symbols:
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
                for signal in comp['signals']:
                    if comp['signals'][signal]['type']=='s_out' and symbols('s_'+signal) in expr.free_symbols:
                        var_=solve(expr, symbols('s_'+signal))
                        solved_expr=var_[0]
                        for ikey,ivar in var_dict_.items():
                            if 'expression' in var_dict_[ikey].keys():
                                solved_expr=solved_expr.subs(symbols(ikey),symbols(var_dict_[ikey]['expression']))
                            else:
                                solved_expr=solved_expr.subs(symbols(ikey),symbols(var_dict_[ikey]['symbol']))
                        comp['constitutive_eqs']+= [[comp['vars']['s_'+signal]['symbol'], ccode(solved_expr), '']]               

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

def update_BG_params(bg_dict, kappa, fName, K, eName, csv_file='params_BG.csv'):
    """
    Update the BG parameters in the bg_dict with the parameters

    Parameters
    ----------
    bg_dict : dict
        The dictionary of the bond graph model
    kappa : numpy.ndarra
        The reaction rate constants
    fName : list
        A list of reaction component Re names
        The order of the names should be the same as the order of the kappa
    K : numpy.ndarray
        The thermodynamic constants
    eName : list
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
        for i in range(len(fName)):
            if abs(kappa[i][0])<1e-2 or abs(kappa[i][0])>1e2:
                kappa_="{:.3e}".format(kappa[i][0])
            else:
                kappa_="{:.3f}".format(kappa[i][0])
            bg_dict[fName[i]]['params']['kappa']['value']=kappa_
            var_name=bg_dict[fName[i]]['params']['kappa']['symbol']
            latex_str=toLatexStr(var_name)
            writer.writerow([latex_str, kappa_])
        for i in range(len(eName)):
            if abs(K[i][0])<1e-2 or abs(K[i][0])>1e2:
                K_="{:.3e}".format(K[i][0])
            else:
                K_="{:.3f}".format(K[i][0])
            bg_dict[eName[i]]['params']['K']['value']=K_
            var_name=bg_dict[eName[i]]['params']['K']['symbol']
            latex_str=toLatexStr(var_name)
            writer.writerow([latex_str, K_])

def EA_BG(bg_dict, bg_ea_components, bg_ea_dict, voi):
    """ Derive the energy analysis and activity equations for the bond graph model

    Parameters
    ----------
    bg_dict : dict
        The dictionary of the bond graph model
    bg_ea_components : dict
        The dictionary of energy analysis for bond graph components
    bg_ea_dict : dict
        The dictionary of energy analysis for the bond graph model
    voi : dict
        The dictionary of the variables of integration
        The default is {
                "description": "Time",
                "units": "second",
                "symbol": "t"
        }

    Returns
    -------
    None
       
    """
    for key, comp in bg_dict.items():
        if 'ports' in comp.keys():
            bg_ea_dict[key]=copy.deepcopy(bg_ea_components[comp['id']]) # bg_ea_components and bg_components have the same id
            bg_ea_dict[key]['constitutive_eqs']=[]
            P_sum=0
            for port in comp['ports']:
                e_str=comp['vars']['e_'+port]['expression'] if 'expression' in comp['vars']['e_'+port].keys() else comp['vars']['e_'+port]['symbol']
                f_str=comp['vars']['f_'+port]['expression'] if 'expression' in comp['vars']['f_'+port].keys() else comp['vars']['f_'+port]['symbol']
                P_symp=sympify(e_str)*sympify(f_str)
                bg_ea_dict[key]['vars']['P_'+port]['symbol']=bg_ea_dict[key]['vars']['P_'+port]['symbol']+ '_' + key
                bg_ea_dict[key]['vars']['P_'+port]['expression']=ccode(P_symp)
                if comp['ports'][port]['direction']=='in':                                     
                    P_sum+=P_symp
                elif comp['ports'][port]['direction']=='out':
                    P_sum+=-P_symp
            bg_ea_dict[key]['vars']['P_sum']['symbol']=bg_ea_dict[key]['vars']['P_sum']['symbol']+ '_' + key
            bg_ea_dict[key]['state_vars']['A']['symbol']=bg_ea_dict[key]['state_vars']['A']['symbol']+ '_' + key
            bg_ea_dict[key]['state_vars']['E']['symbol']=bg_ea_dict[key]['state_vars']['E']['symbol']+ '_' + key
            bg_ea_dict[key]['constitutive_eqs']+=[[bg_ea_dict[key]['vars']['P_sum']['symbol'], ccode(P_sum),'']]
            bg_ea_dict[key]['constitutive_eqs']+=[[bg_ea_dict[key]['state_vars']['A']['symbol'], ccode(Abs(P_sum)),voi['symbol']]]
            bg_ea_dict[key]['constitutive_eqs']+=[[bg_ea_dict[key]['state_vars']['E']['symbol'], ccode(P_sum),voi['symbol']]]      
            
if __name__ == "__main__": 

    file_path='./data/'
    fmatrix='SLC2_f_1.csv'
    rmatrix='SLC2_r_1.csv'
    bg_components_json='BG_components.json'
    bg_components=load_json(bg_components_json)    
    bg_dict={}
    voi={'description': 'Time', 'units': 'second', 'symbol': 't'}   
    direction = 'e2f'
    build_BG_Dict(fmatrix,bg_components,bg_dict,direction,file_path=file_path)
    direction = 'f2e'
    build_BG_Dict(rmatrix,bg_components,bg_dict,direction,file_path=file_path)
    
    update_BG_eqn(bg_dict,voi)
    # dump the bg_dict to a json file, which has the same name and path as the csv file (using Path)
    
    eName, eID, ePort, fName, fID, fPort,N_f=load_matrix(file_path+fmatrix)
    eName, eID, ePort, fName, fID, fPort,N_r=load_matrix(file_path+rmatrix)

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
    update_BG_params(bg_dict, kappa, fName, K, eName, csv_file='params_BG.csv')
   # EA_BG(bg_dict)
    json_file=file_path+'SLC2_BG.json'    
    save_json(bg_dict, json_file)   

    bg_ea_components_json='BG_EA_components.json'
    bg_ea_components=load_json(bg_ea_components_json)
    bg_ea_dict={}
    EA_BG(bg_dict, bg_ea_components, bg_ea_dict,voi)
    json_file=file_path+'SLC2_BG_EA.json'
    save_json(bg_ea_dict, json_file)

        

        
