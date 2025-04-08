# Adapted from https://github.com/BondGraphTools/BondGraphTools 
# under the Apache License http://www.apache.org/licenses/LICENSE-2.0
# Add units to the parameters, variables; 
# Note: Natural logarithm ln=log is used in the constitutive relations

import json

zF = {
        "description": "biochemical-electrical transformer",
        "metamodel": "TF",
        "ports": {
            "0": {
                'orientation': 'in',
                "effort": {
                "description": "chemical potential",
                "units": "J_per_mol",
                "symbol": "mu"  
            }, 
           "flow":{
               "description": "current",
                "units": "fA",
                "symbol": "i"
           },           
                },
            "1": {
                'orientation': 'out',
                "effort": {
                "description": "voltage",
                "units": "volt",
                "symbol": "u" 
            }, 
           "flow":{
                "description": "molar flow",
                "units": "fmol_per_s",
                "symbol": "v"
           },
                },
        },
        "params": {
            "z": {
                "description": "Ratio",
                "value": 1,
                "units": "dimensionless",
                "symbol": "z",
            },
            "F": {
                "description": "Faraday's constant",
                "value": 96485,
                "units": "C_per_mol",
                "symbol": "F",
            },
        },
        "constitutive_relations": [
            "e_1 - z * F * e_0",
            "f_0 - z * F * f_1"
        ]
    }

def BG_components():
    dict_={"description": "biochemical components",
           "domain": "biochemical-electrical",
           
            "physical constants": {
                "F": {
                "description": "Faraday's constant",
                "value": 96485,
                "units": "C_per_mol",
                "symbol": "F",
            },
            },
           "components": {
                "zF": zF,
            }
        }
    return dict_ 

if __name__ == "__main__": 

    json_file = 'EnergyConversion.json'
    comp_dict = BG_components()
    with open(json_file, 'w') as f:
        json.dump(comp_dict, f,indent=4)
