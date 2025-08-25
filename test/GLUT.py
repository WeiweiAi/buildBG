# add the parent directory to the path
import sys
import os
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utilities import load_matrix
from biochemical_params import matrix_constraints, kinetic2BGparams_csvs, BGparams2kinetic_csvs

path_='../data/'

def params_4_state():

    fmatrix=path_+'SLC2_f.csv'
    rmatrix=path_+'SLC2_r.csv'
    eName, eID, ePort, fName, fID, fPort, N_f = load_matrix(fmatrix)
    eName, eID, ePort, fName, fID, fPort, N_r = load_matrix(rmatrix)    
   # V_o=90
   # h=0.726;g=12.1;c=1113;d=90.3;a=500000*V_o;b=a*9.5;f=3000*V_o;e=12.8459*f
    V_i=0.09
    V_o=0.09   
    species=['Ai','Ao', '1','2','3','4']
    volumes={'Ai':V_i,'Ao':V_o}
    constraints=[]
    Ws,K_c,N_c=matrix_constraints(species,volumes,constraints)

    kinetic2BGparams_csvs(K_c,N_c,Ws,eName, fName, N_f, N_r ,kinetic_params_csv=path_+'SLC2_kinetic_params.csv',bg_params_csv=path_+'SLC2_bg_params.csv')
    # convert the BG parameters to the kinetic parameters
    BGparams2kinetic_csvs(Ws,eName,fName, N_f, N_r ,bg_params_csv=path_+'SLC2_bg_params.csv',kinetic_params_csv=path_+'SLC2_kinetic_params_est.csv')



if __name__ == "__main__":
    params_4_state()
