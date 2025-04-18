# Adapted from https://github.com/BondGraphTools/BondGraphTools 
# under the Apache License http://www.apache.org/licenses/LICENSE-2.0
# Add units to the parameters, variables; 
# Adopted Table 1.2 from  Borutzky, Wolfgang. Bond graph modelling of engineering systems. Vol. 103. New York: springer, 2011.
# Note: Natural logarithm ln=log is used in the constitutive relations

import json
Se = {
            "description": "external force source",
            "metamodel": "Se",
            "ports": {
                "0": {
                'orientation': 'out',
                },
           },
            "params": {
                "e": {
                    "description": "force",
                    "units": "J_per_um",
                    "symbol": "F_0",
                    "value": 1
                }
            },
            "constitutive_relations": [
                "e_0 - e"
            ]
        }
Sf = {
            "description": "external velocity source",
            "metamodel": "Sf",
            "ports": {
                "0": {
                'orientation': 'out',
                },
           },
            "params": {
                "f": {
                    "description": "velocity",
                    "units": "um_per_s",
                    "symbol": "v_0",
                    "value": 1
                }
            },
            "constitutive_relations": [
                "f_0 - f"
            ]
        }

R = {
        "description": "linear damper",
        "metamodel": "R",
        "ports": {
             "0": {
                   'orientation': 'in',
                   },
        },
        "params": {
            "R": {
                "description": "coefficient of friction or viscosity",
                "value": 1,
                "units": "J_s_per_um2",
                "symbol": "eta",
            }
        },
        "constitutive_relations": [
            "e_0 - f_0*R"
        ]
    }

C = {
        "description": "linear capacitance",
        "metamodel": "C",
        "ports": {
             "0": {
                  'orientation': 'in',
                   },
        },
        "params": {
            "C": {
                "description": "capacitance",
                "value": 1,
                "units": "um2_per_J",
                "symbol": "C",
            },
        },
        "constitutive_relations": [
            "q_0 - C * e_0",
            "ode(q_0,t) - f_0"
        ]
    }

E = {
        "description": "linear spring",
        "metamodel": "E",
        "ports": {
             "0": {
                  'orientation': 'in',
                   },
        },
        "params": {
            "k": {
                "description": "elastance",
                "value": 1,
                "units": "J_per_um2",
                "symbol": "k",
            },
        },
        "constitutive_relations": [
            "e_0 - k * q_0",
            "ode(q_0,t) - f_0"
        ]
    }

M = {
        "description": "mass",
        "metamodel": "I",
        "ports": {
             "0": {
                  'orientation': 'in',
                   },
        },
        "params": {
            "m": {
                "description": "inertance",
                "value": 1,
                "units": "J_s2_per_um2",
                "symbol": "m",
            },
        },
        "constitutive_relations": [
            "p_0 - m*f_0",
            "ode(p_0,t) - e_0"
        ]
    }

def BG_components():
    dict_={"description": "mechanical components",
           "domain": "mechanical",
           "effort": {
                "description": "force",
                "units": "J_per_um",
                "symbol": "F"  
                }, 
           "flow":{
                "description": "velocity",
                "units": "um_per_s",
                "symbol": "v"
                },
            "generalized displacement":{
                "description": "displacement",
                "units": "um",
                "symbol": "x",
            },
            "generalized momentum":{
                "description": "momentum",
                "units": "J_s_per_um",
                "symbol": "p",         
            },
           "components": {
                "R": R,
                "C": C,
                "E": E,
                "M": M,
                "Se": Se,
                "Sf": Sf,
            }
        }
    return dict_ 

if __name__ == "__main__": 

    json_file = 'BG_mechanical.json'
    comp_dict = BG_components()
    with open(json_file, 'w') as f:
        json.dump(comp_dict, f,indent=4)
