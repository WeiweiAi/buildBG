# Adapted from https://github.com/BondGraphTools/BondGraphTools 
# under the Apache License http://www.apache.org/licenses/LICENSE-2.0
# Add units to the parameters, variables; 
# Adopted Table 1.2 from  Borutzky, Wolfgang. Bond graph modelling of engineering systems. Vol. 103. New York: springer, 2011.
# Note: Natural logarithm ln=log is used in the constitutive relations

import json
Se = {
            "description": "external pressure source",
            "metamodel": "Se",
            "ports": {
                "0": {
                'orientation': 'out',
                },
           },
            "params": {
                "e": {
                    "description": "pressure",
                    "units": "mmHg",
                    "symbol": "P_0",
                    "value": 1
                }
            },
            "constitutive_relations": [
                "e_0 - e"
            ]
        }
Sf = {
            "description": "external flow source",
            "metamodel": "Sf",
            "ports": {
                "0": {
                'orientation': 'out',
                },
           },
            "params": {
                "f": {
                    "description": "volum flow",
                    "units": "mL_per_s",
                    "symbol": "Q_0",
                    "value": 1
                }
            },
            "constitutive_relations": [
                "f_0 - f"
            ]
        }

R = {
        "description": "linear resistance",
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
                "units": "mmHg_s_per_mL",
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
                "units": "mL_per_mmHg_s",
                "symbol": "G",
            }
        },
        "constitutive_relations": [
            "f_0 - e_0*G"
        ]
    }

C = {
        "description": "linear compliance",
        "metamodel": "C",
        "ports": {
             "0": {
                  'orientation': 'in',
                   },
        },
        "params": {
            "C": {
                "description": "volumetric compliance",
                "value": 1,
                "units": "mL_per_mmHg",
                "symbol": "C",
            },
        },
        "constitutive_relations": [
            "q_0 - C * e_0",
            "ode(q_0,t) - f_0"
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
                "description": "volumetric elastance",
                "value": 1,
                "units": "mmHg_per_mL",
                "symbol": "E",
            },
        },
        "constitutive_relations": [
            "e_0 - E * q_0",
            "ode(q_0,t) - f_0"
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
                "description": "inertance",
                "value": 1,
                "units": "mmHg_s2_per_ml",
                "symbol": "L",
            }, 
        },
        "constitutive_relations": [
            "p_0 - L*f_0",
            "ode(p_0, t) - e_0"
        ]
    }

def BG_components():
    dict_={"description": "hydraulic components",
           "domain": "hydraulic",
           "effort": {
                "description": "pressure",
                "units": "mmHg",
                "symbol": "P"  
                }, 
           "flow":{
                "description": "volum flow",
                "units": "mL_per_s",
                "symbol": "Q"
                },
            "generalized displacement":{
                "description": "volume",
                "units": "mL",
                "symbol": "V",
            },
            "generalized momentum":{
                "description": "momentum of a flow tube",
                "units": "mmHg_mL2_per_s3",
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
            }
        }
    return dict_ 

if __name__ == "__main__": 

    json_file = 'BG_hydraulic.json'
    comp_dict = BG_components()
    with open(json_file, 'w') as f:
        json.dump(comp_dict, f,indent=4)
