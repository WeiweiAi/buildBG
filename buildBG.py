from BG_components import e_components_units, biochem_components_units,m_components_units
import pandas as pd
import copy
import numpy as np
import sympy
import csv
from pathlib import Path
import json

defUnit=["ampere","becquerel","candela","celsius","coulomb","dimensionless","farad","gram","gray","henry",
    "hertz","joule","katal","kelvin","kilogram","liter","litre","lumen","lux","meter","metre","mole",
    "newton","ohm","pascal","radian","second","siemens","sievert","steradian","tesla","volt","watt","weber"]
params_common=['R','T','F']

def load_matrix(fmatrix,rmatrix):
    """
    Load stoichiometric matrices from csv files

    Parameters
    ----------
    fmatrix : str
        The file path of the forward stoichiometric matrix
    rmatrix : str
        The file path of the reverse stoichiometric matrix

    Returns
    -------
    CompName : list
        A list of component names
    CompType : list
        A list of component types
    ReName : list
        A list of reaction names
    ReType : list
        A list of reaction types
    N_f : numpy.ndarray
        The forward stoichiometric matrix
    N_r : numpy.ndarray
        The reverse stoichiometric matrix
    """
    # * * ReType ReType
    # * * ReName ReName
    # CompType CompName 0 1 
    # CompType CompName 1 0
    startR=2
    startC=2
    N_f = []
    N_r = []
    CompName=[]
    CompType=[]
    with open(fmatrix,'r') as f:
        reader = csv.reader(f,delimiter=',')
        line_count = 0
        for row in reader:
            if line_count ==startR-2:
                ReType=row[startC:]
                line_count += 1
            elif line_count ==startR-1:
                ReName=[s for s in row[startC:]]
                line_count += 1
            else:
                N_f.append(row[startC:])
                CompName.append(row[startC-1])
                CompType.append(row[startC-2])
        f.close()
    with open(rmatrix,'r') as f:
        reader = csv.reader(f,delimiter=',')
        line_count = 0
        for row in reader:
            if line_count <startR:                
                line_count += 1            
            else:
                N_r.append(row[startC:])
        f.close()
    
    return CompName, CompType, ReName, ReType, np.array(N_f).astype(int), np.array(N_r).astype(int)

def kinetic2BGparams(N_f,N_r,kf,kr,K_c,N_c,Ws):
    """
    Convert kinetic parameters to BG parameters

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
        The volume vector,
        a column vector with the size of number of reactions nr + number of species ns,
        the first nr elements are 1 for reactions, the last ns elements are the volume of species

    Returns
    -------
    kappa : numpy.ndarray
        The reaction rate constants
    K : numpy.ndarray
        The thermodynamic constants
    K_eq : numpy.ndarray
        The equilibrium constants
    diff : float
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
    W=np.vstack([np.ones((num_cols,1)),Ws])   
    # convert kinetic parameters to BG parameters
    lambdaW= np.exp(np.matmul(np.linalg.pinv(M),np.log(k)))
    lambda_ = np.divide(lambdaW,W)
    kappa=lambda_[:num_cols]
    K = lambda_[num_cols:]
    
    # check if the solution is valid
    N_rref, _ = sympy.Matrix(N).rref()
    zero_est = None
    R_mat = np.array(sympy. nsimplify(sympy.Matrix(N), rational=True).nullspace())
    if R_mat.size>0:
        R_mat = np.transpose(np.array(R_mat).astype(np.float64))[0]
        zero_est = np.matmul(R_mat.T,K_eq)
    # Check that there is a detailed balance constraint
    Z = sympy.nsimplify(sympy.Matrix(N_b), rational=True).nullspace() #rational_nullspace(M, 2)
    if Z:
        Z = np.transpose(np.array(Z).astype(np.float64))[0]
        zero_est = np.matmul(Z.T,np.log(K_contraints))

    k_est = np.exp(np.matmul(M,np.log(lambdaW)))
    diff = np.sum(np.abs(np.divide(k_est - k,k)))

    return kappa, K, K_eq, diff, zero_est

def buildBG(fmatrix,rmatrix,file_path='./'):
    e_components=e_components_units()['components']
    biochem_components=biochem_components_units()['components']
    m_components=m_components_units()['components']
    CompName,CompType,ReName,ReType,N_f,N_r=load_matrix(file_path+fmatrix,file_path+rmatrix)
    compNames=CompName+ReName
    compTypes=CompType+ReType
    n_zeros=len(CompName)
    n_ones=len(ReName)
    # Use the CompType and ReType to look up the corresponding components in the e_components and biochem_components
    # Get the parameters for components and declare the parameters in the format: "var varName: varUnits {init: 1, pub: out}}"
    # Save the parameters in a dictionary with the key as the variable name and the value as the cellml code
    comp_dict={}    
    for i in range(len(compNames)):
        compType=compTypes[i]
        compName=compNames[i]
        compIndex=compName # compName is the key for the dictionary, make it unique; the ReName and CompName may have overlap;
        if compType in e_components.keys():
            comp_dict[compIndex]=copy.deepcopy(e_components[compType])
            comp_dict[compIndex]['type']=compType
            # Instantiate the parameters, variables and constitutive relations for the component
            for param in comp_dict[compIndex]['params'].keys():
                if param not in params_common:
                    comp_dict[compIndex]['params'][param]['symbol']=comp_dict[compIndex]['params'][param]['symbol']+ '_' + compName
            for var in comp_dict[compIndex]['vars'].keys():
                comp_dict[compIndex]['vars'][var]['symbol']=comp_dict[compIndex]['vars'][var]['symbol']+ '_' + compName
            if 'state_vars' in e_components[compType].keys():                  
                for state_var in e_components[compType]['state_vars'].keys():
                    comp_dict[compIndex]['state_vars'][state_var]['symbol']=comp_dict[compIndex]['state_vars'][state_var]['symbol']+ '_' + compName                  
        elif compType in biochem_components.keys():
            comp_dict[compIndex]=copy.deepcopy(biochem_components[compType])
            comp_dict[compIndex]['type']=compType
            # Instantiate the parameters, variables and constitutive relations for the component
            for param in comp_dict[compIndex]['params'].keys():
                if param not in params_common:
                    comp_dict[compIndex]['params'][param]['symbol']=comp_dict[compIndex]['params'][param]['symbol']+ '_' + compName
            for var in comp_dict[compIndex]['vars'].keys():
                comp_dict[compIndex]['vars'][var]['symbol']=comp_dict[compIndex]['vars'][var]['symbol']+ '_' + compName
            if 'state_vars' in biochem_components[compType].keys():                    
                for state_var in biochem_components[compType]['state_vars'].keys():
                    comp_dict[compIndex]['state_vars'][state_var]['symbol']=comp_dict[compIndex]['state_vars'][state_var]['symbol']+ '_' + compName
        elif compType in m_components.keys():
            comp_dict[compIndex]=copy.deepcopy(m_components[compType])
            comp_dict[compIndex]['type']=compType
            # Instantiate the parameters, variables and constitutive relations for the component
            for param in comp_dict[compIndex]['params'].keys():
                if param not in params_common:
                    comp_dict[compIndex]['params'][param]['symbol']=comp_dict[compIndex]['params'][param]['symbol']+ '_' + compName
            for var in comp_dict[compIndex]['vars'].keys():
                comp_dict[compIndex]['vars'][var]['symbol']=comp_dict[compIndex]['vars'][var]['symbol']+ '_' + compName
            if 'state_vars' in m_components[compType].keys():                    
                for state_var in m_components[compType]['state_vars'].keys():
                    comp_dict[compIndex]['state_vars'][state_var]['symbol']=comp_dict[compIndex]['state_vars'][state_var]['symbol']+ '_' + compName
        else:
            print('The component type is not found in the e_components or biochem_components or m_components')

    for j in range(len(ReName)):
        # The e_0 of the R component is the sum of the e_0 of each C component in the column of N_f[i,:]
        # The e_1 of the R component is the sum of the e_0 of each C component in the column of N_r[i,:]
        reIndex=ReName[j]
        for i in range(len(CompName)):
            compIndex=CompName[i]
            # 0 node to 1 node on port 0; 
            if N_f[i,j]!=0:
                comp_dict[reIndex]['ports']['0']['in']+=[compIndex+f':{N_f[i,j]}']
            # 0 node to 1 node on port 1; 
            if N_r[i,j]!=0:
                comp_dict[reIndex]['ports']['1']['in']+=[compIndex+f':{-N_r[i,j]}']
            # 1 node to 0 node on port 0; 
            if  'f_0' in comp_dict[compIndex]['vars'].keys() and comp_dict[compIndex]['vars']['f_0']['IOType']=='in':
                if N_f[i,j]!=0:
                    comp_dict[compIndex]['ports']['0']['in']+=[reIndex+f':{-N_f[i,j]}']

                if N_r[i,j]!=0:
                    comp_dict[compIndex]['ports']['0']['in']+=[reIndex+f':{N_r[i,j]}']
            # 1 node to 0 node on port 1;
            if 'f_1' in comp_dict[compIndex]['vars'].keys() and comp_dict[compIndex]['vars']['f_1']['IOType']=='in':
                if N_f[i,j]!=0:
                    comp_dict[compIndex]['ports']['1']['in']+=[reIndex+f':{-N_f[i,j]}']

                if N_r[i,j]!=0:
                    comp_dict[compIndex]['ports']['1']['in']+=[reIndex+f':{N_r[i,j]}']

    update_eqn(comp_dict)
    return comp_dict

def update_params(comp_dict,n_zeros, kappa, K, q_init_all, csv='params_BG.csv'):
    # assume that the kappa and K are in the same order as the components in the comp_dict
    # Create a pd frame with the columns: Parameter, Value, and Unit
    csv_pd=pd.DataFrame(columns=['Parameter','Value','Unit'])
    if K.size>0:
        j=0
        for i in range(len(K)):
            compIndex=list(comp_dict)[i]
            if abs(K[i][0])<1e-2 or abs(K[i][0])>1e2:
                K_="{:.3e}".format(K[i][0])
            else:
                K_="{:.3f}".format(K[i][0])
            comp_dict[compIndex]['params']['K']['value']=K_
            if 'q_init' in comp_dict[compIndex]['params'].keys() and comp_dict[compIndex]['type']=='Ce':
                comp_dict[compIndex]['params']['q_init']['value']=q_init_all[j][0]
                j+=1
            str_var=comp_dict[compIndex]['params']['K']['symbol']
            # split the string with the underscore and get the first element
            sub_str_1=str_var.split('_')[0]
            if len(str_var.split('_'))>1:
                sub_str_2=str_var.split('_')[1]
                latex_str='$'+sub_str_1+'_{'+sub_str_2+'}$'
            else:
                sub_str_2=''
                latex_str='$'+sub_str_1+'$'

            units=comp_dict[compIndex]['params']['K']['units']
            # split the string with the underscore and 'per'
            unit_latex_str=''
            if 'per' in units:
                unit_latex_str+=str.join('.',units.split('per')[0].split('_'))
                units_per=units.split('per')[1]
                for j in range(len(units_per.split('_'))-1):
                    unit_latex_str+=units_per.split('_')[j+1]+'$^{-1}$'
            else:
                unit_latex_str+=str.join('.',units.split('_'))

            csv_pd.loc[i]=[latex_str,K_,unit_latex_str]
    if kappa.size>0:
        for i in range(len(kappa)):
            compIndex=list(comp_dict)[i+n_zeros]
            if abs(kappa[i][0])<1e-2 or abs(kappa[i][0])>1e2:
                kappa_="{:.3e}".format(kappa[i][0])
            else:
                kappa_="{:.3f}".format(kappa[i][0])
            comp_dict[compIndex]['params']['kappa']['value']=kappa_
            str_var=comp_dict[compIndex]['params']['kappa']['symbol']
            # split the string with the underscore and get the first element
            sub_str_1=str_var.split('_')[0]
            if len(str_var.split('_'))>1:
                sub_str_2=str_var.split('_')[1]
                latex_str='$\\'+sub_str_1+'_{'+sub_str_2+'}$'
            else:
                sub_str_2=''
                latex_str='$\\'+sub_str_1+'$'
            units=comp_dict[compIndex]['params']['kappa']['units']
            # split the string with the underscore and 'per'
            unit_latex_str=''
            if 'per' in units:
                unit_latex_str+=str.join('.',units.split('per')[0].split('_'))
                units_per=units.split('per')[1]
                for j in range(len(units_per.split('_'))-1):
                    unit_latex_str+=units_per.split('_')[j+1]+'$^{-1}$'
            else:
                unit_latex_str+=str.join('.',units.split('_'))
           
            csv_pd.loc[i+len(K)]=[latex_str,kappa_,unit_latex_str]

    csv_pd.to_csv(csv,index=False)
    # dump the comp_dict to a json file, which has the same name and path as the csv file (using Path)
    json_file=Path(csv).with_suffix('.json')    
    with open(json_file, 'w') as f:
        json.dump(comp_dict, f,indent=4)    

def update_eqn(comp_dict):
    def get_flow_outputs(sub_comps):
        flow_outputs=[]
        for comp in sub_comps:
            comp_name=comp.split(':')[0]
            stochoimetry=comp.split(':')[1]
            try:
                num_stochoimetry=int(stochoimetry)
            except:
                try:
                    num_stochoimetry=float(stochoimetry)
                except:
                    raise ValueError('The stochoimetry is not an integer or a float')

            for var in comp_dict[comp_name]['vars']:
                if 'Flow' in comp_dict[comp_name]['vars'][var]['description'] and comp_dict[comp_name]['vars'][var]['IOType']=='out':
                    if stochoimetry=='1':
                        flow_outputs+=[f"+{comp_dict[comp_name]['vars'][var]['symbol']}"]
                    elif stochoimetry=='-1':
                        flow_outputs+=[f"-{comp_dict[comp_name]['vars'][var]['symbol']}"]
                    elif num_stochoimetry>0:
                        flow_outputs+=[f"+{num_stochoimetry}{{dimensionless}}*{comp_dict[comp_name]['vars'][var]['symbol']}"]
                    else:
                        flow_outputs+=[f"-{-num_stochoimetry}{{dimensionless}}*{comp_dict[comp_name]['vars'][var]['symbol']}"]
                        
        return flow_outputs 
    def get_efforts_outputs(sub_comps):
        efforts_outputs=[]
        for comp in sub_comps:
            comp_name=comp.split(':')[0]
            stochoimetry=comp.split(':')[1]
            try:
                num_stochoimetry=int(stochoimetry)
            except:
                try:
                    num_stochoimetry=float(stochoimetry)
                except:
                    raise ValueError('The stochoimetry is not an integer or a float')
            for var in comp_dict[comp_name]['vars']:
                if 'Potential' in comp_dict[comp_name]['vars'][var]['description'] and comp_dict[comp_name]['vars'][var]['IOType']=='out':
                    if stochoimetry=='1':
                        efforts_outputs+=[f"+{comp_dict[comp_name]['vars'][var]['symbol']}"]
                    elif stochoimetry=='-1':
                        efforts_outputs+=[f"-{comp_dict[comp_name]['vars'][var]['symbol']}"]
                    elif num_stochoimetry>0:
                        efforts_outputs+=[f"+{num_stochoimetry}{{dimensionless}}*{comp_dict[comp_name]['vars'][var]['symbol']}"]
                    else:
                        efforts_outputs+=[f"-{-num_stochoimetry}{{dimensionless}}*{comp_dict[comp_name]['vars'][var]['symbol']}"]
        return efforts_outputs
    
    for key, comp in comp_dict.items():
        if comp['type']=='Ce':
            comp['constitutive_relations']=[f"{comp['vars']['e_0']['symbol']} = R*T*ln({comp['params']['K']['symbol']}*{comp['state_vars']['q_0']['symbol']})",
            f"ode({comp['state_vars']['q_0']['symbol']},t) = {comp['vars']['f_0']['symbol']}"]
            # Get all in flows
            in_flows_=comp['ports']['0']['in']
            in_flows_vars=get_flow_outputs(in_flows_)
            sum_in_flows=' '.join(in_flows_vars)
            sum_in_flows= sum_in_flows[1:] if sum_in_flows[0]=='+' else sum_in_flows
            comp['conservation_relations']=[
                f"{comp['vars']['f_0']['symbol']} = {sum_in_flows}",
                ]               
        elif comp['type']=='ch_Se':
            comp['constitutive_relations']=[f"{comp['vars']['e_0']['symbol']} = R*T*ln({comp['params']['K']['symbol']}*{comp['params']['q_init']['symbol']})"]            
        elif comp['type']=='Re':
            comp['constitutive_relations']=[
            f"{comp['vars']['f_0']['symbol']} = {comp['params']['kappa']['symbol']}*(exp({comp['vars']['e_0']['symbol']}/(R*T)) - exp({comp['vars']['e_1']['symbol']}/(R*T)))"
          ]
             # Get all in effors
            in_efforts_port0_=comp['ports']['0']['in']
            in_efforts_port0_vars=get_efforts_outputs(in_efforts_port0_)
            sum_in_efforts_port0=' '.join(in_efforts_port0_vars)
            sum_in_efforts_port0= sum_in_efforts_port0[1:] if sum_in_efforts_port0[0]=='+' else sum_in_efforts_port0

            in_efforts_port1_=comp['ports']['1']['in']
            in_efforts_port1_vars=get_efforts_outputs(in_efforts_port1_)
            sum_in_efforts_port1=' '.join(in_efforts_port1_vars)
            sum_out_efforts_port1 = sum_in_efforts_port1.replace('+','plus').replace('-','+').replace('plus','-')
            sum_out_efforts_port1= sum_out_efforts_port1[1:] if sum_out_efforts_port1[0]=='+' else sum_out_efforts_port1
            
            comp['conservation_relations']=[
                f"{comp['vars']['e_0']['symbol']} = {sum_in_efforts_port0}",
                f"{comp['vars']['e_1']['symbol']} = {sum_out_efforts_port1}",
            ]
        elif comp['type']=='C' or comp['type']=='m_C':
            comp['constitutive_relations']=[
            f" {comp['vars']['e_0']['symbol']}={comp['state_vars']['q_0']['symbol']}/{comp['params']['C']['symbol']}",
            f"ode({comp['state_vars']['q_0']['symbol']},t) = {comp['vars']['f_0']['symbol']}"
          ]
             # Get all in flows
            in_flows_=comp['ports']['0']['in']
            in_flows_vars=get_flow_outputs(in_flows_)
            sum_in_flows=' '.join(in_flows_vars)
            sum_in_flows= sum_in_flows[1:] if sum_in_flows[0]=='+' else sum_in_flows
            comp['conservation_relations']=[
                f"{comp['vars']['f_0']['symbol']} = {sum_in_flows}",
                ]

        elif comp['type']=='R' or comp['type']=='m_R':
            comp['constitutive_relations']=[
                f"{comp['vars']['f_0']['symbol']}={comp['vars']['e_0']['symbol']}/{comp['params']['r']['symbol']}",
            ]
            # Get all in effors
            in_efforts_=comp['ports']['0']['in']
            in_efforts_vars=get_efforts_outputs(in_efforts_)
            sum_in_efforts=' '.join(in_efforts_vars)
            sum_in_efforts= sum_in_efforts[1:] if sum_in_efforts[0]=='+' else sum_in_efforts
            comp['conservation_relations']=[
                f"{comp['vars']['e_0']['symbol']} = {sum_in_efforts}",
                ]             
        elif comp['type']=='Se' or comp['type']=='e_Se' or comp['type']=='m_Se':
            comp['constitutive_relations']=[
                f"{comp['vars']['e_0']['symbol']}={comp['params']['e']['symbol']}"
            ]
        elif comp['type']=='Sf' or comp['type']=='e_Sf' or comp['type']=='m_Sf':
            comp['constitutive_relations']=[
                f"{comp['vars']['f_0']['symbol']}={comp['params']['f']['symbol']}"
            ]
        elif comp['type']=='TF':
            comp['constitutive_relations']=[
                f"{comp['vars']['e_1']['symbol']}={comp['params']['r']['symbol']}*{comp['vars']['e_0']['symbol']}",
                f"{comp['vars']['f_0']['symbol']}={comp['params']['r']['symbol']}*{comp['vars']['f_1']['symbol']}"
            ]
            # Get all in flows on port 1
            in_flows_=comp['ports']['1']['in']
            in_flows_vars=get_flow_outputs(in_flows_)
            sum_in_flows=' '.join(in_flows_vars)
            sum_in_flows= sum_in_flows[1:] if sum_in_flows[0]=='+' else sum_in_flows
            comp['conservation_relations']=[
                f"{comp['vars']['f_1']['symbol']} = {sum_in_flows}",
            ]                      
            # Get all in effors on port 0
            in_efforts_=comp['ports']['0']['in']
            in_efforts_vars=get_efforts_outputs(in_efforts_)
            sum_in_efforts=' '.join(in_efforts_vars)
            sum_out_efforts=sum_in_efforts.replace('+','plus').replace('-','+').replace('plus','-')
            sum_out_efforts= sum_out_efforts[1:] if sum_out_efforts[0]=='+' else sum_out_efforts
            comp['conservation_relations']+=[
                f"{comp['vars']['e_0']['symbol']} = {sum_out_efforts}"
                ]
        elif comp['type']=='zF':
            comp['constitutive_relations']=[
                f"{comp['vars']['e_1']['symbol']}={comp['params']['r']['symbol']}*F*{comp['vars']['e_0']['symbol']}",
                f"{comp['vars']['f_0']['symbol']}={comp['params']['r']['symbol']}*F*{comp['vars']['f_1']['symbol']}"
            ]
             # Get all in flows on port 1
            in_flows_=comp['ports']['1']['in']
            in_flows_vars=get_flow_outputs(in_flows_)
            sum_in_flows=' '.join(in_flows_vars)
            sum_in_flows= sum_in_flows[1:] if sum_in_flows[0]=='+' else sum_in_flows
            comp['conservation_relations']=[
                f"{comp['vars']['f_1']['symbol']} = {sum_in_flows}",
            ]                      
            # Get all in effors on port 0
            in_efforts_=comp['ports']['0']['in']
            in_efforts_vars=get_efforts_outputs(in_efforts_)
            sum_in_efforts=' '.join(in_efforts_vars)
            sum_out_efforts=sum_in_efforts.replace('+','plus').replace('-','+').replace('plus','-')
            sum_out_efforts= sum_out_efforts[1:] if sum_out_efforts[0]=='+' else sum_out_efforts
            comp['conservation_relations']+=[
                f"{comp['vars']['e_0']['symbol']} = {sum_out_efforts}"
                ]

        elif comp['type']=='GY':
            comp['constitutive_relations']=[
                f"{comp['vars']['e_1']['symbol']}={comp['params']['r']['symbol']}*{comp['vars']['f_0']['symbol']}",
                f"{comp['vars']['e_0']['symbol']}={comp['params']['r']['symbol']}*{comp['vars']['f_1']['symbol']}"
            ]
    return comp_dict

def to_cellmlV1_params(comp_dict, model_name='params_BG',model_file='params_BG.txt',file_path='./'):
    indent=' '*4
    cellml_code=f'def model {model_name} as\n'
    cellml_code+=indent+'def import using "./units.cellml" for\n'
    param_units=set()
    params_common_set=set()
    # Get all the units used in the parameters
    for comp in comp_dict:
        for param in comp_dict[comp]['params']:
            if param in params_common:
                params_common_set.add(param)
            if comp_dict[comp]['params'][param]['units'] not in defUnit:
                param_units.add(comp_dict[comp]['params'][param]['units'])
    for unit in param_units:
        cellml_code+=indent*2+f"unit {unit} using unit {unit};\n"
    cellml_code+=indent+'enddef;\n'

    cellml_code+=indent+f'def comp {model_name} as\n'
    if 'R' in params_common_set:
        cellml_code+=indent*2+f"var R: J_per_K_mol"+f"{{ init: 8.31, pub: out}};\n"
    if 'T' in params_common_set:
        cellml_code+=indent*2+f"var T: kelvin"+f"{{ init: 293, pub: out}};\n"
    if 'F' in params_common_set:
        cellml_code+=indent*2+f'var F: C_per_mol'+f"{{ init: 96485, pub: out}};\n"
    for comp in comp_dict:
        for param in comp_dict[comp]['params']:
            if param not in ['R','T','F']:
                if (comp_dict[comp]['type']=="ch_Se" and param=='q_init') or (comp_dict[comp]['type']=='e_Se' and param =='e') or (comp_dict[comp]['type']=="m_Se" and param=='e'): 
                    cellml_code+=indent*2+f"var {comp_dict[comp]['params'][param]['symbol']}: {comp_dict[comp]['params'][param]['units']}" + f"{{ pub: out}};\n"
                else:
                    cellml_code+=indent*2+f"var {comp_dict[comp]['params'][param]['symbol']}: {comp_dict[comp]['params'][param]['units']}" + f"{{ init: {comp_dict[comp]['params'][param]['value']}, pub: out}};\n"
                    
        
    cellml_code+=indent+'enddef;\n'
    cellml_code+='enddef;\n'
    # Save the model to a file
    with open(file_path+model_file, 'w') as f:
        f.write(cellml_code)
    return cellml_code

def to_cellmlV1_models(comp_dict, model_name='BG',model_file='BG.txt',params_file='params_BG.cellml',file_path='./'):
    indent=' '*4
    # Apart from parameters, add vairables, state variables, constitutive relations and conservation relations
    param_units=set()
    cellml_code_p=''
    cellml_code_vars=''
    cellml_code_state_vars=''
    cellml_code_constitutive_relations=''
    cellml_code_conservation_relations=''
    params_common_set=set()
    for comp in comp_dict:
        for param in comp_dict[comp]['params']:
            if param in params_common:
                params_common_set.add(param)
            if comp_dict[comp]['params'][param]['units'] not in defUnit:
                param_units.add(comp_dict[comp]['params'][param]['units'])
        for var in comp_dict[comp]['vars']:
            if comp_dict[comp]['vars'][var]['units'] not in defUnit:
                param_units.add(comp_dict[comp]['vars'][var]['units'])
        if 'state_vars' in comp_dict[comp].keys():
            for state_var in comp_dict[comp]['state_vars']:
                param_units.add(comp_dict[comp]['state_vars'][state_var]['units'])
    
    cellml_code=f'def model {model_name} as\n'
    cellml_code+=indent+'def import using "./units.cellml" for\n'
    for unit in param_units:
        cellml_code+=indent*2+f"unit {unit} using unit {unit};\n"
    cellml_code+=indent+'enddef;\n'
    
    cellml_code+=indent+f'def import using "{params_file}" for\n'
    params_model_name=params_file.split('.')[0]
    cellml_code+=indent*2+f'comp {params_model_name} using comp {params_model_name};\n'
    cellml_code+=indent+'enddef;\n'
    
    cellml_code+=indent+f'def comp {model_name} as\n'
    cellml_code+=indent*2+f"var t: second;\n"
    if 'R' in params_common_set:
        cellml_code+=indent*2+f"var R: J_per_K_mol"+f"{{ pub: in}};\n"
    if 'T' in params_common_set:
        cellml_code+=indent*2+f"var T: kelvin"+f"{{ pub: in}};\n" 
    if 'F' in params_common_set:
        cellml_code+=indent*2+f'var F: C_per_mol'+f"{{ pub: in}};\n"
    for comp in comp_dict:
        for param in comp_dict[comp]['params']:
            if param not in params_common:
                cellml_code_p+=indent*2+f"var {comp_dict[comp]['params'][param]['symbol']}: {comp_dict[comp]['params'][param]['units']}" + f"{{ pub: in}};\n"
        for var in comp_dict[comp]['vars']:
            cellml_code_vars+=indent*2+f"var {comp_dict[comp]['vars'][var]['symbol']}: {comp_dict[comp]['vars'][var]['units']};\n" 
        if 'state_vars' in comp_dict[comp].keys():
            for state_var in comp_dict[comp]['state_vars']:
                q_init=comp_dict[comp]['state_vars'][state_var]['value']
                cellml_code_state_vars+=indent*2+f"var {comp_dict[comp]['state_vars'][state_var]['symbol']}: {comp_dict[comp]['state_vars'][state_var]['units']}" + f"{{ init: {comp_dict[comp]['params'][q_init]['symbol']}}};\n"
        for enq in comp_dict[comp]['constitutive_relations']:
            cellml_code_constitutive_relations+=indent*2+enq+';\n'
        if 'conservation_relations' in comp_dict[comp].keys():
            for enq in comp_dict[comp]['conservation_relations']:
                cellml_code_conservation_relations+=indent*2+enq+';\n'
    cellml_code+=cellml_code_p
    cellml_code+=cellml_code_vars
    cellml_code+=cellml_code_state_vars
    cellml_code+=cellml_code_constitutive_relations
    cellml_code+=cellml_code_conservation_relations
    cellml_code+=indent+f"enddef;\n"
    cellml_code+=indent+f'def map between {params_model_name} and {model_name} for\n'
    if 'R' in params_common_set:
        cellml_code+=indent*2+'vars R and R;\n'
    if 'T' in params_common_set:
        cellml_code+=indent*2+'vars T and T;\n'
    if 'F' in params_common_set:
        cellml_code+=indent*2+'vars F and F;\n'
    for comp in comp_dict:
        for param in comp_dict[comp]['params']:
            if param not in params_common:
                cellml_code+=indent*2+f"vars {comp_dict[comp]['params'][param]['symbol']} and {comp_dict[comp]['params'][param]['symbol']};\n"
    cellml_code+=indent+'enddef;\n'

    cellml_code+='enddef;\n'
    # Save the model to a file
    with open(file_path+model_file, 'w') as f:
        f.write(cellml_code)
    return cellml_code

            
if __name__ == "__main__": 
    
    file_path='./data/'
    fmatrix='SLC2_f.csv'
    rmatrix='SLC2_r.csv'
    comp_dict=buildBG(fmatrix,rmatrix,file_path)
    CompName,CompType,ReName,ReType,N_f,N_r=load_matrix(file_path+fmatrix,file_path+rmatrix)

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

    kappa, K, K_eq, diff,  zero_est= kinetic2BGparams(N_f,N_r,kf,kr,K_c,N_c,Ws)
    n_zeros=len(CompName)

    AVO=6.022e23
    q_tot=6e10 # 5e10 molecules per cell
    q_init = q_tot/AVO/6*1e15 # fmol
    q_init_all=np.array([[q_init]*4]).transpose() 

    csv_file=file_path+'SLC2_BG.csv'
    update_params(comp_dict,n_zeros, kappa, K, q_init_all, csv_file)
    to_cellmlV1_params(comp_dict, model_name='params_BG',model_file='params_BG.txt',file_path=file_path)
    to_cellmlV1_models(comp_dict, model_name='GLUT2_BG',model_file='GLUT2_BG.txt',params_file='params_BG.cellml',file_path=file_path)

        

        
