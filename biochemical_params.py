from math import exp
from sympy import nsimplify,Matrix
import numpy as np
import csv
import os
from utilities import load_matrix
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
    kappa : 1d numpy.ndarray
        The reaction rate constants
    K : 1d numpy.ndarray
        The thermodynamic constants
    Ws : 1d list, the size is the number of species
           
    Returns
    -------
    kf : 1d list
        The forward rate constants, the same order as the reactions in N_f  
    kr : 1d list
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
    # update the BG parameters for the biochemical reactions
    eName, eID, ePort, fName, fID, fPort,N_f=load_matrix(fmatrix)
    eName, eID, ePort, fName, fID, fPort,N_r=load_matrix(rmatrix)
    """
    var{environment.K_1} K_1: per_fmol {init: 101619537.2009};
    var{environment.K_10} K_10: per_fmol {init: 20459.5509};
    var{environment.K_11} K_11: per_fmol {init: 121.4456};
    var{environment.K_12} K_12: per_fmol {init: 3.1436};
    var{environment.K_13} K_13: per_fmol {init: 0.32549};
    var{environment.K_14} K_14: per_fmol {init: 156.3283};
    var{environment.K_15} K_15: per_fmol {init: 1977546.8577};
    var{environment.K_2} K_2: per_fmol {init: 63209.8623};
    var{environment.K_3} K_3: per_fmol {init: 157.2724};
    var{environment.K_4} K_4: per_fmol {init: 14.0748};
    var{environment.K_5} K_5: per_fmol {init: 5.0384};
    var{environment.K_6} K_6: per_fmol {init: 92.6964};
    var{environment.K_7} K_7: per_fmol {init: 4854.5924};
    var{environment.K_8} K_8: per_fmol {init: 15260.9786};
    var{environment.K_9} K_9: per_fmol {init: 13787022.8009};
    var{environment.K_H} K_H: per_fmol {init: 0.04565};
    var{environment.K_Ke} K_Ke: per_fmol {init: 0.009236};
    var{environment.K_Ki} K_Ki: per_fmol {init: 0.0012595};
    var{environment.K_MgADP} K_MgADP: per_fmol {init: 7.976e-05};
    var{environment.K_MgATP} K_MgATP: per_fmol {init: 2.3715};
    var{environment.K_Nae} K_Nae: per_fmol {init: 0.0061242};
    var{environment.K_Nai} K_Nai: per_fmol {init: 0.00083514};
    var{environment.K_P} K_P: per_fmol {init: 0.04565};
    var{environment.kappa_1} kappa_1: fmol_per_sec {init: 330.5462};
    var{environment.kappa_10} kappa_10: fmol_per_sec {init: 259461.6507};
    var{environment.kappa_11} kappa_11: fmol_per_sec {init: 172042.3334};
    var{environment.kappa_12} kappa_12: fmol_per_sec {init: 6646440.3909};
    var{environment.kappa_13} kappa_13: fmol_per_sec {init: 597.4136};
    var{environment.kappa_14} kappa_14: fmol_per_sec {init: 70.9823};
    var{environment.kappa_15} kappa_15: fmol_per_sec {init: 0.015489};
    var{environment.kappa_2} kappa_2: fmol_per_sec {init: 132850.9145};
    var{environment.kappa_3} kappa_3: fmol_per_sec {init: 200356.0223};
    var{environment.kappa_4} kappa_4: fmol_per_sec {init: 2238785.3951};
    var{environment.kappa_5} kappa_5: fmol_per_sec {init: 10787.9052};
    var{environment.kappa_6} kappa_6: fmol_per_sec {init: 15.3533};
    var{environment.kappa_7} kappa_7: fmol_per_sec {init: 2.3822};
    var{environment.kappa_8} kappa_8: fmol_per_sec {init: 2.2855};
    var{environment.kappa_9} kappa_9: fmol_per_sec {init: 1540.1349};
    """
    kappa_1= 330.5462
    kappa_2= 132850.9145
    kappa_3= 200356.0223
    kappa_4= 2238785.3951
    kappa_5= 10787.9052
    kappa_6= 15.3533
    kappa_7= 2.3822
    kappa_8= 2.2855
    kappa_9= 1540.1349
    kappa_10= 259461.6507
    kappa_11= 172042.3334
    kappa_12= 6646440.3909
    kappa_13= 597.4136
    kappa_14= 70.9823
    kappa_15= 0.015489
    kappa=[kappa_1,kappa_2,kappa_3,kappa_4,kappa_5,kappa_6,kappa_7,kappa_8,kappa_9,kappa_10,kappa_11,kappa_12,kappa_13,kappa_14,kappa_15]
    K_1= 101619537.2009
    K_2= 63209.8623
    K_3= 157.2724
    K_4= 14.0748
    K_5= 5.0384
    K_6= 92.6964
    K_7= 4854.5924
    K_8= 15260.9786
    K_9= 13787022.8009
    K_10= 20459.5509
    K_11= 121.4456
    K_12= 3.1436
    K_13= 0.32549
    K_14= 156.3283
    K_15= 1977546.8577
    K_Ki=0.0012595
    K_Ke=0.009236
    K_Nai=0.00083514
    K_Nae=0.0061242
    K_MgATP=2.3715
    K_MgADP=7.976e-05
    K_P=0.04565
    K_H=0.04565
    K=[K_1,K_2,K_3,K_4,K_5,K_6,K_7,K_8,K_9,K_10,K_11,K_12,K_13,K_14,K_15,K_Ki,K_Ke,K_Nai,K_Nae,K_MgATP,K_MgADP,K_P,K_H]
    W_i=38
    W_e=5.182
    Ws=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,W_i,W_e,W_i,W_e,W_i, W_i, W_i, W_i]
    k_f,k_r=BGparams2kinetic(N_f,N_r,kappa,K,Ws)
    a=K_MgATP/(K_MgADP*K_H*K_P)
    print('a:',a)
    k_dict={'k_f':k_f,'k_r':k_r}
    k_df=pd.DataFrame(data=k_dict)
    k_df.to_csv(file_path+'kinetic_params.csv', )
    
    k_f_name=['k_f_'+str(i+1) for i in range(len(k_f))]
    k_r_name=['k_r_'+str(i+1) for i in range(len(k_r))]
    # Plot the k_f and k_r values in a bar chart to see the normalized values; no seaborn, only matplotlib
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 6))
    bar_width = 0.35
    index = np.arange(len(k_f))  # Create tick locations for the reactions

    # Plot the bars for k_f and k_r
    ax.bar(index, k_f, bar_width, label='k_f', color='b')
    ax.bar(index + bar_width, k_r, bar_width, label='k_r', color='r')

    # Set the x-axis labels
    ax.set_xlabel('Reaction rate constants')
    ax.set_ylabel('Values')
    ax.set_title('Kinetic Parameters')

    # Set the tick locations and labels
    tick_locations = np.concatenate([index, index + bar_width])  # Combine tick locations for k_f and k_r
    tick_labels = k_f_name + k_r_name  # Combine labels for k_f and k_r
    ax.set_xticks(index + bar_width / 2)  # Center the ticks between the bars
    ax.set_xticklabels(k_f_name, rotation=45)  # Use only `k_f_name` for the x-axis labels

    # Add a legend and save the plot
    ax.legend()
    plt.tight_layout()
    plt.savefig(file_path + 'kinetic_params.png')
    
    # convert the kinetic parameters to BG parameters
    deltaG_ATP=11.9e3 #J/mol
    R=8.314 #J/(K*mol)
    T=310 #K
    K_c=[exp(-deltaG_ATP/(R*T))*10**6, 1,1,K_Ki*W_i,K_MgATP*W_i,K_MgADP*W_i,K_P*W_i,K_H*W_i] # the equilibrium constant for ATP hydrolysis
    N_c=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,-1],
         [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0],
         [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0],
         [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
         [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
         [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
         [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
         [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]]
    kappa_new, K_new, K_eq, diff_,  zero_est,k_est= kinetic2BGparams(N_f,N_r,k_f,k_r,K_c,N_c,Ws)
    print('zero_est:',zero_est)
    print('diff:',diff_)
    k_est=k_est.flatten()
    num_k=len(k_est)-len(K_c)
    print('num_k:',num_k)
    k_dict_new={'k_f':k_est[:int(num_k/2)],'k_r':k_est[int(num_k/2):30]}
    k_df_new=pd.DataFrame(data=k_dict_new)
    k_df_new.to_csv(file_path+'kinetic_params_new.csv', )
    param_val=np.concatenate((kappa_new,K_new))
    param_name=fName+eName
    param_units=['$fmol/s$']*len(fName)+['$fmol^{-1}$']*len(eName)
    write_params_csv(param_name,param_val,param_units,csv_file=file_path+'params_BG.csv')
    

    
        
    