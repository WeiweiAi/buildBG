# Adapted from https://github.com/BondGraphTools/BondGraphTools 
# under the Apache License http://www.apache.org/licenses/LICENSE-2.0
# Add units to the parameters, variables; 
# Adopted Table 1.2 from  Borutzky, Wolfgang. Bond graph modelling of engineering systems. Vol. 103. New York: springer, 2011.
# Note: Natural logarithm ln=log is used in the constitutive relations

import json
R = {
        "description": "linear resistor",
        "metamodel": "R",
        "ports": {
            "0": {
                'orientation': 'in',
                },
        },
        "params": {
            "R": {
                "description": "resistance",
                "value": 1,
                "units": "pOhm",
                "symbol": "R",
            },
        },        
        "constitutive_relations": [
            "e_0 - f_0*R",
        ]
    }

G = {
        "description": "linear conductance",
        "metamodel": "G",
        "ports": {
            "0": {
                'orientation': 'in',
                },
        },
        "params": {
            "G": {
                "description": "conductance",
                "value": 1,
                "units": "tS",
                "symbol": "G",
            },
        },        
        "constitutive_relations": [
            "f_0 - e_0*G"
        ]
    }

C = {
        "description": "linear capacitor",
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
                "units": "fF",
                "symbol": "C",
            },
        },
        "constitutive_relations": [
            "q_0 - C * e_0",
            "ode(q_0,t) - f_0",
        ]
    }

E = {
        "description": "linear elastance",
        "metamodel": "E",
        "ports": {
            "0": {
                'orientation': 'in',
                },
        },
        "params": {
            "E": {
                "description": "elastance",
                "value": 1,
                "units": "per_fF",
                "symbol": "E",
            },
        },
        "constitutive_relations": [
            "e_0 - E * q_0",
            "ode(q_0,t) - f_0",
        ]
    }

I = {
        "description": "linear inductor",
        "metamodel": "I",
        "ports": {
            "0": {
                'orientation': 'in',
                },
        },
        "params": {
            "L": {
                "description": "inductance",
                "value": 1,
                "units": "pH",
                "symbol": "L",
            },
        },
        "constitutive_relations": [
            "p_0 - L*f_0",
            "ode(p_0,t) - e_0",
        ]
    }

Se = {
        "description": "external voltage source",
        "metamodel": "Se",
         "ports": {
            "0": {
                'orientation': 'out',
                },
        },
        "params": {
            "e": {
                "description": "voltage",
                "units": "volt",
                "symbol": "u_0",
                "value": 1
            }
        },
        "constitutive_relations": [
            "e_0 - e",
        ]
    }


Sf = {
        "description": "external current source",
        "metamodel": "Sf",
        "ports": {
            "0": {
                'orientation': 'out',
                },
        },
        "params": {
            "f": {
                "description": "current",
                "units": "fA",
                "symbol": "i_0",
                "value": 1
            }
        },
        "constitutive_relations": [
            "f_0 - f",
        ]
    }

TF = {
        "description": "linear transformer",
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
            "n": {
                "description": "ratio",
                "value": 1,
                "units": "dimensionless",
                "symbol": "n",
            }
        },

        "constitutive_relations": [
            "e_1 - n * e_0",
            "f_0 - n * f_1",
        ]
    }

def BG_components():
    dict_={"description": "electrical components",
           "domain": "electrical",
           "effort": {
                "description": "voltage",
                "units": "volt",
                "symbol": "u" 
            }, 
           "flow":{
               "description": "current",
                "units": "fA",
                "symbol": "i"
           },
            "generalized displacement":{
                "description": "charge",
                "units": "fC",
                "symbol": "q",
            },
            "generalized momentum":{
                "description": "magnetic flux linkage",
                "units": "volt_s",
                "symbol": "p",         
            },
           "components": {
                "R": R,
                "G": G,
                "C": C,
                "E": E,
                "I": I,
                "Se": Se,
                "Sf": Sf,
                "TF": TF,
            }
        }
    return dict_ 

if __name__ == "__main__": 

    json_file = 'BG_electrical.json'
    comp_dict = BG_components()
    with open(json_file, 'w') as f:
        json.dump(comp_dict, f,indent=4)
