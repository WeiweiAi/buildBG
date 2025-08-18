from sympy import  nsimplify,Matrix
import numpy as np
import csv
import os
from utilities import load_matrix, load_matrix_domain
import pandas as pd

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

def BGparams2kinetic(N_f,N_r,kappa,K,Ws):
    """
    Convert the BG parameters to the kinetic parameters

    Parameters
    ----------
    N_f : 2d numpy.ndarray
        The forward stoichiometry matrix
    N_r : 2d numpy.ndarray
        The reverse stoichiometry matrix
    kappa : 1d list
        The reaction rate constants
    K : 1d list
        The thermodynamic constants
    Ws : 1d list, the size is the number of species
           
    Returns
    -------
    kf : 1d numpy.ndarray
        The forward rate constants, the same order as the reactions in N_f  
    kr : 1d numpy.ndarray
        The reverse rate constants, the same order as the reactions in N_r
    
    """

    N_fT=np.transpose(N_f)
    N_rT=np.transpose(N_r)
    N = N_r - N_f
    num_cols = N_f.shape[1] # number of reactions, the same as the number of columns in N_f
    num_rows = N_f.shape[0] # number of species, the same as the number of rows in N_f
    I=np.identity(num_cols)
    M=np.block([
        [I, N_fT],
        [I, N_rT]
    ])
    W_s=np.array([Ws]).transpose()
    W=np.vstack([np.ones((num_cols,1)),W_s])
    kappa_ = np.array([kappa]).transpose()
    K_ = np.array([K]).transpose()
    lambda_= np.vstack([kappa_,K_])
    lambdaW = np.multiply(lambda_,W)
    k = np.exp(np.matmul(M,np.log(lambdaW)))
    k_f = k[:num_cols]
    k_r = k[num_cols:]
    return k_f[:,0],k_r[:,0]

def kinetic2BGparams_csvs(Kc,Nc,Ws,eName, fName, N_f, N_r, kinetic_params_csv='kinetic_params.csv',bg_params_csv='bg_params.csv'):
    """
    Convert the BG parameters to the kinetic parameters

    Parameters
    ----------
    Kc : 1d list, the constraints, can be empty
        The constraints, can be empty
    Nc : 2d list, can be empty, 
        the length of the Nc is the number of Kc,
        the length of the Nc[0] is the number of the species
    Ws : 1d list, the size is the number of species
    eName : 1d list, the size is the number of species
    fName : 1d list, the size is the number of reactions
    N_f: 2d numpy.ndarray
        The forward stoichiometry matrix
    N_r: 2d numpy.ndarray
        The reverse stoichiometry matrix
    kinetic_params_csv : str, optional
        The file path of the kinetic parameters csv file
        The default is 'kinetic_params.csv'.
    bg_params_csv : str, optional
        The file path of the bg parameters csv file
        The default is 'bg_params.csv'.

    returns
    -------
    None

    side effect
    ------------
    Save the bg parameters to a csv file
    Save the estimated kinetic parameters to a new csv file            
    """
    
    # check if the number of species is the same as the number of rows in the forward stoichiometry matrix
    if len(Ws)!=N_f.shape[0]:
        raise ValueError('The number of species in the volume csv file is not the same as the number of rows in the forward stoichiometry matrix')
    pdf_k= pd.read_csv(kinetic_params_csv,index_col=0)
    k_f=pdf_k['k_f'].to_numpy().flatten()
    k_r=pdf_k['k_r'].to_numpy().flatten()
    # the number of reactions is the same as the number of columns in the forward stoichiometry matrix
    if len(k_f)!=N_f.shape[1]:
        raise ValueError('The number of reactions in the kinetic parameters csv file is not the same as the number of columns in the forward stoichiometry matrix')
    if len(k_r)!=N_r.shape[1]:
        raise ValueError('The number of reactions in the kinetic parameters csv file is not the same as the number of columns in the reverse stoichiometry matrix')
    kappa, K, K_eq, diff_, zero_est,k_est= kinetic2BGparams(N_f,N_r,k_f,k_r,Kc,Nc,Ws)
    print('zero_est:',zero_est)
    print('diff:',diff_)
    bg_dict = {}
    # create a dictionary for the bg parameters
    for i in range(len(kappa)):
        key=f'kappa_{fName[i]}'
        bg_dict[key] = kappa[i]
    for i in range(len(K)):
        key=f'K_{eName[i]}'
        bg_dict[key] = K[i]
    # save the dictionary to a csv file
    pd.DataFrame.from_dict(bg_dict, orient='index').to_csv(bg_params_csv, header=False)
    # save the recalculated k_est parameters to a csv file
    k_est=k_est.flatten()
    num_k=len(k_est)-len(Kc)
    k_dict_new={'k_f':k_est[:int(num_k/2)],'k_r':k_est[int(num_k/2):num_k]}
    k_df_new=pd.DataFrame(data=k_dict_new)
    # index starts from 1
    k_df_new.index = [i+1 for i in range(len(k_f))]
    new_csv_file = kinetic_params_csv.split('.csv')[0] + '_new.csv'
    k_df_new.to_csv(new_csv_file, index_label='Reaction')
     
def BGparams2kinetic_csvs(Ws,eName,fName, N_f, N_r,bg_params_csv='bg_params.csv',kinetic_params_csv='kinetic_params.csv'):
    """
    Convert the BG parameters to the kinetic parameters

    Parameters
    ----------
    Ws : 1d list, the size is the number of species
    eName : 1d list, the size is the number of species
    fName : 1d list, the size is the number of reactions
    N_f : 2d numpy.ndarray
        The forward stoichiometry matrix
    N_r : 2d numpy.ndarray
        The reverse stoichiometry matrix
    bg_params_csv : str, optional
        The file path of the bg parameters csv file
        The default is 'bg_params.csv'.
    kinetic_params_csv : str, optional
        The file path of the kinetic parameters csv file
        The default is 'kinetic_params.csv'.

    Returns
    -------
    None

    side effect
    ------------
    Save the kinetic parameters to a csv file

    """
    # check if the number of species is the same as the number of rows in the forward stoichiometry matrix
    if len(Ws)!=N_f.shape[0]:
        raise ValueError('The number of species in the volume csv file is not the same as the number of rows in the forward stoichiometry matrix')
    bg_dict = pd.read_csv(bg_params_csv, header=None).set_index(0).T.to_dict('records')[0]
    # get K and kappa from the bg_dict
    kappa = [bg_dict[f'kappa_{f}'] for f in fName]
    K = [bg_dict[f'K_{e}'] for e in eName]
    # check if the number of reactions is the same as the number of columns in the forward stoichiometry matrix
    if len(kappa)!=N_f.shape[1]:
        raise ValueError('The number of reactions in the bg parameters csv file is not the same as the number of columns in the forward stoichiometry matrix')
    if len(K)!=N_r.shape[0]:
        raise ValueError('The number of reactions in the bg parameters csv file is not the same as the number of rows in the reverse stoichiometry matrix')
    # convert the bg parameters to the kinetic parameters
    k_f,k_r=BGparams2kinetic(N_f,N_r,kappa,K,Ws)
    # save the kinetic parameters to a csv file
    k_dict={'k_f':k_f,'k_r':k_r}
    k_df=pd.DataFrame(data=k_dict)
    # index starts from 1
    k_df.index = [i+1 for i in range(len(k_f))]
    k_df.to_csv(kinetic_params_csv, index_label='Reaction')    

def matrix_constraints(species,volumes, constraints=None):
    """
    Create the constraints matrix and volume vector for the BG parameters

    parameters
    ----------
    species : list
        The list of species
    volumes : dict
        The dictionary of volumes and the keys are the species names
    constraints : dict, optional
        The dictionary of constraints,[{'num':[species],'denom':[species],'value':1}]
    """
    # default volumes
    Ws = np.ones(len(species))
    # update volumes
    for i, s in enumerate(species):
        if s in volumes.keys():
            Ws[i] = volumes[s]
    K_c = []
    N_c = []
    if constraints:
        cols=len(constraints)
        for i in range(cols):
            num = [] if 'num' not in constraints[i].keys() else constraints[i]['num']
            denom = [] if 'denom' not in constraints[i].keys() else constraints[i]['denom']
            value = constraints[i]['value']
            K_c.append(value)
            nc=[0]*len(species)
            # for the species in the numerator, get their indices and set the corresponding list element 1;
            for s in num:
                nc[species.index(s)] = 1
            for s in denom:
                nc[species.index(s)] = -1
            N_c.append(nc)

    return Ws,K_c, N_c

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

if __name__ == "__main__": 
    current_dir = os.path.dirname(os.path.abspath(__file__))
    file_path=os.path.join(current_dir, './data/')
    fmatrix=file_path+'NKE15state_f_chem.csv'
    rmatrix=file_path+'NKE15state_r_chem.csv'
    
    K_Ki=0.0012595
    K_Ke=0.009236
    K_Nai=0.00083514
    K_Nae=0.0061242
    K_MgATP=2.3715
    K_MgADP=7.976e-05
    K_P=0.04565
    K_H=0.04565
    W_i=38
    W_e=5.182
    K_1=101619537.2009
    # convert the kinetic parameters to BG parameters
    deltaG_ATP=11.9e3 #J/mol
    R=8.314 #J/(K*mol)
    T=310 #K
    species=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','Ki','Ke','Nai','Nae','MgATP','MgADP','P','H']
    volumes={'Ki':W_i,'Ke':W_e,'Nai':W_i,'Nae':W_e,'MgATP':W_i,'MgADP':W_i,'P':W_i,'H':W_i}
    constraints=[{'num':['MgATP'],'denom':['MgADP','P','H'],'value':np.exp(-deltaG_ATP/(R*T))*10**6},
                 {'num':['Ki'],'value':K_Ki*W_i},
                 {'num':['MgADP'],'value':K_MgADP*W_i},
                 {'num':['P'],'value':K_P*W_i},
                 {'num':['H'],'value':K_H*W_i},
                 {'num':['MgATP'],'value':K_MgATP*W_i},
                 {'num':['1'],'value':K_1}]

    Ws,K_c,N_c=matrix_constraints(species,volumes,constraints)
    eDomain, eName, eID, ePort, fDomain, fName, fID, fPort, N_f = load_matrix_domain(fmatrix)
    eDomain, eName, eID, ePort, fDomain, fName, fID, fPort, N_r = load_matrix_domain(rmatrix)

    kinetic2BGparams_csvs(K_c,N_c,Ws,eName, fName, N_f, N_r ,kinetic_params_csv=file_path+'kinetic_params.csv',bg_params_csv=file_path+'bg_params.csv')
    # convert the BG parameters to the kinetic parameters
    BGparams2kinetic_csvs(Ws,eName,fName, N_f, N_r ,bg_params_csv=file_path+'bg_params.csv',kinetic_params_csv=file_path+'kinetic_params_est.csv')

    
        
    