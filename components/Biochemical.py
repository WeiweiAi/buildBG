# Adapted from https://github.com/BondGraphTools/BondGraphTools 
# under the Apache License http://www.apache.org/licenses/LICENSE-2.0
# Add units to the parameters, variables; 
# Note: Natural logarithm ln=log is used in the constitutive relations

import json
Ce = {
        "description":"chemical species",
        "metamodel":"C",
        "ports": {
            "0": {
                'orientation': 'in',
                },
        },
        "params":{
          "K":{
              "description": "thermodynamic constant; exp(mu_0/RT)/V",
              "units": "per_fmol",
              "symbol": "K",
              "value": 1
          },
        
        },    
        "constitutive_relations":[
          "e_0 - R*T*log(K*q_0)",
          "ode(q_0,t) - f_0"
        ]
     }

Re ={
        "description": "biochemical reaction",
        "metamodel":"R",
        "ports": {
            "0": {
                'orientation': 'in',
                },
            "1": {
                'orientation': 'out',
                },
        },
        "params":{
          "kappa":{
              "description":"reaction rate",
              "units": "fmol_per_s",
              "value": 1,
              "symbol": "kappa"
          },
        },
       
        "constitutive_relations":[
          "f_0 - kappa*(exp(e_0/R/T) - exp(e_1/R/T))",
          "f_0 - f_1"
        ]
    }

zF = {
        "description": "biochemical-electrical transformer",
        "metamodel": "TF",
        "ports": {
            "0": {
                'orientation': 'in',
                },
            "1": {
                'orientation': 'out',
                },
        },
        "params": {
            "z": {
                "description": "Ratio",
                "value": 1,
                "units": "dimensionless",
                "symbol": "z",
            },
        },
        "constitutive_relations": [
            "e_1 - z * F * e_0",
            "f_0 - z * F * f_1"
        ]
    }

Re_GHK ={
        "description": "voltage modulated biochemical reaction-Goldman-Hodgkin-Katz (GHK) ion channel",
        "metamodel":"MR",
        "ports": {
            "0": {
                'orientation': 'in',
                },
            "1": {
                'orientation': 'out',
                },
        },
        "signals": {
             "m": {
                "description": "modulation port",
                'orientation': 'in',
                "units": "J_per_mol",
                "symbol": "Am",
                },
        },
        "params":{
          "kappa":{
              "description":"reaction rate",
              "units": "fmol_per_s",
              "value": 1,
              "symbol": "kappa"
          },
        },
        "constitutive_relations":[
          "f_0- Piecewise(( kappa*(exp(e_0/R/T) - exp(e_1/R/T)),s_m==0),(kappa*(s_m/R/T/(exp(s_m/R/T)))*((exp(e_0/R/T) - exp(e_1/R/T))),s_m!=0))",
          "f_0-f_1"
        ]
    }

Se = {
          "description":"chemostat",
          "metamodel":"Se",
          "ports": {
            "0": {
                'orientation': 'in',
                },
         },
          "params":{
            "K":{
                "description": "thermodynamic constant; exp(mu_0/RT)/V",
                "units": "per_fmol",
                "symbol": "K",
                "value": 1
            },
            "q_0":{
                "description":"molar quantity",
                "units": "fmol",
                "symbol": "q_0",
                "value": 1,
            },
          },
          "constitutive_relations":[
            "e_0 - R*T*log(K*q_0)",
          ]
        }

Sf = {
        "description": "flow source",
        "metamodel": "Sf",
        "ports": {
            "0": {
                'orientation': 'out',
                },
        },
        "params": {
            "f": {
                "description": "molar flow",
                "units": "fmol_per_s",
                "symbol": "v_0",
                "value": 1
            }
        },
        "constitutive_relations": [
            "f_0 - f",
        ]
    }

def BG_components():
    dict_={"description": "biochemical components",
           "domain": "biochemical",
           "effort": {
                "description": "chemical potential",
                "units": "J_per_mol",
                "symbol": "mu"  
            }, 
           "flow":{
                "description": "molar flow",
                "units": "fmol_per_s",
                "symbol": "v"
           },
           "generalized displacement":{
                "description":"molar quantity",
                "units": "fmol",
                "symbol": "q",
            },
            "physical constants": {
                "R":{
                "description":"universal gas constant",
                "units": "J_per_K_mol",
                "symbol": "R",
                "value": 8.31,
                },
                "F": {
                "description": "Faraday's constant",
                "value": 96485,
                "units": "C_per_mol",
                "symbol": "F",
            },
            },
            "thermodynamic parameters":{
                "T":{
                "description": "temperature",
                "units": "kelvin",
                "symbol": "T",
                "value": 293,
                },
            },
           "components": {
                "Ce": Ce,
                "Re": Re,
                "Re_GHK": Re_GHK,
                "zF": zF,
                "Se": Se,
                "Sf": Sf
            }
        }
    return dict_ 

if __name__ == "__main__": 

    json_file = 'BG_biochemical.json'
    comp_dict = BG_components()
    with open(json_file, 'w') as f:
        json.dump(comp_dict, f,indent=4)
