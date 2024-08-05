# Adapted from https://github.com/BondGraphTools/BondGraphTools 
# under the Apache License http://www.apache.org/licenses/LICENSE-2.0
# Add units to the parameters, variables; Add power and activity variables; Add signals; 
# Define individual components
# Note: Natural logarithm ln=log is used in the constitutive relations
import json
R = {
        "description": "Generalised Linear Resistor",
        "domain": "Electrical",
        "id": "R",
        "metamodel": "R",
        "ports": {
            "0": {'in':[],'out':[],'type': ' ','direction': 'in'},
        },
        "signals": {
        },
        "params": {
            "r": {
                "description": "Resistance",
                "value": 1,
                "units": "pOhm",
                "symbol": "R",
                "range": [
                    0,
                    "inf"
                ]
            }
        },
        "vars":{
            "e_0": {
                "description": "Generalised Potential",
                "units": "volt",
                "symbol": "V",                
            },
            "f_0": {
                "description": "Generalised Flow",
                "units": "fA",
                "symbol": "I",
            },
        },
        "constitutive_relations": [
            "e_0 - f_0*r",

        ]
    }
R_EA = {'componentID': 'R',
        "vars":{
            "P_0": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P",
            },
            "P_sum": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P_sum",
            }
        },
        "state_vars": {
            "A": {
                "description": "Generalised activity",
                "units": "fJ",
                "symbol": "A",
                "value": 0,
            },
            "E": {
                "description": "Generalised energy",
                "units": "fJ",
                "symbol": "E",
                "value": 0,
            }
        },
    }
C = {
        "description": "Generalised Linear Capacitor",
        "domain": "Electrical",
        "id": "C",
        "metamodel": "C",
        "ports": {
            "0": {'in':[],'out':[],'type': ' ','direction': 'in'}
        },
        "signals": {
        },
        "params": {
            "C": {
                "description": "Capacitance",
                "value": 1,
                "units": "fF",
                "symbol": "C",
                "range": [
                    0,
                    "inf"
                ]
            },
            "q_init": {
                "description": "Initial Generalised Position",
                "units": "fC",
                "symbol": "q_init",     
                "value": 1,
            },
        },
        "vars":{
            "e_0": {
                "description": "Generalised Potential",
                "units": "volt",
                "symbol": "V",
            },
            "f_0": {
                "description": "Generalised Flow",
                "units": "fA",
                "symbol": "I",
            },
        },
        
        "state_vars": {
            "q_0": {
                "description": "Generalised position",
                "units": "fC",
                "symbol": "q",
                "value": "q_init"
            },
        },
        "constitutive_relations": [
            "q_0 - C * e_0",
            "ode(q_0,t) - f_0",
        ]
    }
C_EA = {'componentID': 'C',
        "vars":{
            "P_0": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P",
            },
            "P_sum": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P_sum",
            }
        },
        "state_vars": {
            "A": {
                "description": "Generalised activity",
                "units": "fJ",
                "symbol": "A",
                "value": 0,
            },
            "E": {
                "description": "Generalised energy",
                "units": "fJ",
                "symbol": "E",
                "value": 0,
            }
        },
    }
I = {
        "description": "Generalised Linear Inductor",
        "domain": "Electrical",
        "id": "I",
        "metamodel": "I",
        "ports": {
            "0": {'in':[],'out':[],'type': ' ','direction': 'in'}
        },
        "signals": {
        },
        "vars":{
            "e_0": {
                "description": "Generalised Potential",
                "units": "volt",
                "symbol": "V",
            },
            "f_0": {
                "description": "Generalised Flow",
                "units": "fA",
                "symbol": "I",
            },
        },
        "state_vars": {
            "p_0": {
                "description": "Generalised momentum",
                "units": "volt_s",
                "symbol": "p",
                "value": "p_init"
            }
        },
        "params": {
            "L": {
                "description": "Inductance",
                "value": 1,
                "units": "pH",
                "symbol": "L",
                "range": [
                    0,
                    "inf"
                ]
            },
            "p_init": {
                "description": "Initial generalised momentum",
                "units": "volt_s",
                "symbol": "p_init",
                "value": 1
            }
        },
        "constitutive_relations": [
            "p_0 - L*f_0",
            "ode(p_0,t) - e_0",
        ]
    }
I_EA = {'componentID': 'I',
        "vars":{
            "P_0": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P",
            },
            "P_sum": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P_sum",
            }
        },
        "state_vars": {
            "A": {
                "description": "Generalised activity",
                "units": "fJ",
                "symbol": "A",
                "value": 0,
            },
            "E": {
                "description": "Generalised energy",
                "units": "fJ",
                "symbol": "E",
                "value": 0,
            }
        },
    }
Se = {
        "description": "Effort Source",
        "domain": "Electrical",
        "id": "Se",
        "metamodel": "SS",
        "ports": {
            "0": {'in':[],'out':[],'type': ' ','direction': 'out'}
        },
        "signals": {
        },
        "params": {
            "e": {
                "description": "Generalised Potential",
                "units": "volt",
                "symbol": "V0",
                "value": 0
            }
        },
        "vars":{
            "e_0": {
                "description": "Generalised Potential",
                "units": "volt",
                "symbol": "V",
            },
            "f_0": {
                "description": "Generalised Flow",
                "units": "fA",
                "symbol": "I",
            },
            
        },
        "constitutive_relations": [
            "e_0 - e",
        ]
    }

Se_EA = {'componentID': 'Se',
        "vars":{
            "P_0": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P",
            },
            "P_sum": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P_sum",
            }
        },
        "state_vars": {
            "A": {
                "description": "Generalised activity",
                "units": "fJ",
                "symbol": "A",
                "value": 0,
            },
            "E": {
                "description": "Generalised energy",
                "units": "fJ",
                "symbol": "E",
                "value": 0,
            }
        },
}

Sf = {
        "description": "Flow Source",
        "domain": "Electrical",
        "id": "Sf",
        "metamodel": "SS",
        "ports": {
            "0": {'in':[],'out':[],'type': ' ','direction': 'out'}
        },
        "signals": {
        },
        "params": {
            "f": {
                "description": "Generalised Flow",
                "units": "fA",
                "symbol": "I0",
                "value": 0
            }
        },
        "vars":{
            "e_0": {
                "description": "Generalised Potential",
                "units": "volt",
                "symbol": "V",
            },
            "f_0": {
                "description": "Generalised Flow",
                "units": "fA",
                "symbol": "I",
            },
        },
        "constitutive_relations": [
            "f_0 - f",
        ]
    }

Sf_EA = {'componentID': 'Sf',
        "vars":{
            "P_0": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P",
            },
            "P_sum": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P_sum",
            }
        },
        "state_vars": {
            "A": {
                "description": "Generalised activity",
                "units": "fJ",
                "symbol": "A",
                "value": 0,
            },
            "E": {
                "description": "Generalised energy",
                "units": "fJ",
                "symbol": "E",
                "value": 0,
            }
        },
}   

TF = {
        "description": "Linear Transformer",
        "domain": "Electrical",
        "id": "TF",
        "metamodel": "TF",
        "ports": {
            "0": {
                "description": "Primary",
                'in':[],'out':[],'type': ' ','direction': 'in'
            },
            "1": {
                "description": "Secondary",
                'in':[],'out':[],'type': ' ','direction': 'out'
            }
        },
        "signals": {
        },
        "vars":{
            "e_0": {
                "description": "Generalised Potential",
                "units": "volt",
                "symbol": "V_0",
            },
            "f_0": {
                "description": "Generalised Flow",
                "units": "fA",
                "symbol": "I_0",
            },
            "e_1": {
                "description": "Generalised Potential",
                "units": "volt",
                "symbol": "V_1",
            },
            "f_1": {
                "description": "Generalised Flow",
                "units": "fA",
                "symbol": "I_1",
            },
        },
        "params": {
            "r": {
                "description": "Ratio",
                "value": 1,
                "units": "dimensionless",
                "symbol": "n",
            }
        },

        "constitutive_relations": [
            "e_1 - r * e_0",
            "f_0 - r * f_1",
        ]
    }

TF_EA = {'componentID': 'TF',
        "vars":{
            "P_0": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P_0",
            },
            "P_1": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P_1",
            },
            "P_sum": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P_sum",
            }
        },
        "state_vars": {
            "A": {
                "description": "Generalised activity",
                "units": "fJ",
                "symbol": "A",
                "value": 0,
            },
            "E": {
                "description": "Generalised energy",
                "units": "fJ",
                "symbol": "E",
                "value": 0,
            }
        },
}

GY = {
        "description": "Linear Gyrator",
        "domain": "Electrical-Mechanical",
        "id": "GY",
        "metamodel": "GY",
        "ports": {
            "0": {
                "description": "Primary",
                'in':[],'out':[],'type': ' ','direction': 'in'
            },
            "1": {
                "description": "Secondary",
                'in':[],'out':[],'type': ' ','direction': 'out'
            }
        },
        "signals": {
        },
        "vars":{
            "e_0": {
                "description": "Generalised Potential",
            },
            "f_0": {
                "description": "Generalised Flow",
            },
            "e_1": {
                "description": "Generalised Potential",
            },
            "f_1": {
                "description": "Generalised Flow",
            },

        },
        "params": {
            "r": {
                "description": "Ratio",
                "value": 1
            }
        },
        "constitutive_relations": [
            "e_1 - r*f_0",
            "e_0 - r*f_1"
        ]
    }

Ce = {
        "description":"Concentration of Chemical Species",
        "domain": "biochemical",
        "id":"Ce",
        "metamodel":"C",
        "ports":{
          "0":{
              'in':[],'out':[],
              'type': ' ','direction': 'in'
          },
        },
        "signals": {
        },
        "params":{
          "K":{
              "description": "Biochemical Constant; exp(mu_0/RT)/V",
              "units": "per_fmol",
              "symbol": "K",
              "value": 1
          },
          "R":{
              "description":"Universal Gas Constant",
              "units": "J_per_K_mol",
              "symbol": "R",
              "value": 8.31,
          },
          "T":{
              "description": "Temperature",
              "units": "kelvin",
              "symbol": "T",
              "value": 293,
          },
          "q_init":{
              "description": "Initial Molar Quantity",
              "units": "fmol",
              "symbol": "q_init",
              "value": 1,
          }
        },    
        "state_vars":{ 
          "q_0":{
              "description":"Molar Quantity",
              "units": "fmol",
              "symbol": "q",
              "value": "q_init",
          },
        },
        "vars":{
          "e_0": {
              "description": "Generalised Potential",
              "units": "J_per_mol",
              "symbol": "mu",
          },
          "f_0": {
              "description": "Generalised Flow ",
              "units": "fmol_per_s",
              "symbol": "v",
          }      
        },
        "constitutive_relations":[
          "e_0 - R*T*log(K*q_0)",
          "ode(q_0,t) - f_0"
        ]
     }

Ce_EA = {'componentID': 'Ce',
        "vars":{
            "P_0": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P",
            },
            "P_sum": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P_sum",
            }
        },
        "state_vars": {
            "A": {
                "description": "Generalised activity",
                "units": "fJ",
                "symbol": "A",
                "value": 0,
            },
            "E": {
                "description": "Generalised energy",
                "units": "fJ",
                "symbol": "E",
                "value": 0,
            }
        },
}

Re ={
        "description": "Biochemical Reaction",
        "domain": "biochemical",
        "id": "Re",
        "metamodel":"R",
        "ports":{
          "0":{'in':[],'out':[],'type': ' ','direction': 'in'},
          "1":{'in':[],'out':[],'type': ' ','direction': 'out'}
        },
        "signals": {
        },
        "params":{
          "kappa":{
              "description":"Reaction Rate",
              "units": "fmol_per_s",
              "value": 1,
              "symbol": "kappa"
          },
          "R":{
              "description":"Universal Gas Constant",
              "value": 8.31,
              "units": "J_per_K_mol",
              "symbol": "R"
          },
          "T":{
              "description": "Temperature",
              "value": 293,
              "units": "kelvin",
              "symbol": "T"
          }
        },
        "vars":{
              "e_0": {
                  "description": "Generalised Potential",
                  "units": "J_per_mol",
                  "symbol": "Af",
              },
              "f_0": {
                  "description": "Generalised Flow",
                  "units": "fmol_per_s",
                  "symbol": "v",
              },
              "e_1": {
                  "description": "Generalised Potential",
                  "units": "J_per_mol",
                  "symbol": "Ar",
              },
              "f_1": {
                  "description": "Generalised Flow",
                  "units": "fmol_per_s",
                  "symbol": "v",
                  "expression": "v" # exclude it from the variable list
              },
        },
        "constitutive_relations":[
          "f_0 - kappa*(exp(e_0/R/T) - exp(e_1/R/T))"
        ]
    }

Re_EA = {'componentID': 'Re',
        "vars":{
             "P_0": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P_0",
              },
              "P_1": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P_1",
             },
            "P_sum": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P_sum",
            }
        },
        "state_vars": {
            "A": {
                "description": "Generalised activity",
                "units": "fJ",
                "symbol": "A",
                "value": 0,
            },
            "E": {
                "description": "Generalised energy",
                "units": "fJ",
                "symbol": "E",
                "value": 0,
            }
        },
}

zF = {
        "description": "Biochemical Transformer",
        "domain": "biochemical",
        "id": "zF",
        "metamodel": "TF",
        "ports": {
            "0": {
                "description": "Primary",
                'in':[],'out':[],'type': ' ','direction': 'in'
            },
            "1": {
                "description": "Secondary",
                'in':[],'out':[],'type': ' ','direction': 'out'
            }
        },
        "signals": {
        },
        "vars":{
            "e_0": {
                "description": "Generalised Potential",
                "units": "volt",
                "symbol": "V",
            },
            "f_0": {
                "description": "Generalised Flow",
                "units": "fA",
                "symbol": "I",
            },
            "e_1": {
                "description": "Generalised Potential",
                "units": "J_per_mol",
                "symbol": "mu",
            },
            "f_1": {
                "description": "Generalised Flow",
                "units": "fmol_per_s",
                "symbol": "v",
            },
        },
        "params": {
            "r": {
                "description": "Ratio",
                "value": 1,
                "units": "dimensionless",
                "symbol": "z",
            },
            "F": {
                "description": "Faraday's Constant",
                "value": 96485,
                "units": "C_per_mol",
                "symbol": "F",
            }
        },
        "constitutive_relations": [
            "e_1 - r * F* e_0",
            "f_0 - r * F * f_1"
        ]
    }

zF_EA = {'componentID': 'zF',
         "vars":{
                "P_0": {
                    "description": "Generalised power",
                    "units": "fJ_per_s",
                    "symbol": "P_0",
                },
                "P_1": {
                    "description": "Generalised power",
                    "units": "fJ_per_s",
                    "symbol": "P_1",
                },
            "P_sum": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P_sum",
            }
            },
        "state_vars": {
            "A": {
                "description": "Generalised activity",
                "units": "fJ",
                "symbol": "A",
                "value": 0,
            },
            "E": {
                "description": "Generalised energy",
                "units": "fJ",
                "symbol": "E",
                "value": 0,
            }
        },
}

ch_Se = {
          "description":"Concentration of Chemical Species",
          "domain": "biochemical",
          "id":"ch_Se",
          "metamodel":"SS",
          "ports":{
            "0":{
                'in':[],'out':[],
                'type': ' ','direction': 'out'
            },
          },
          "signals": {
        },
          "params":{
            "K":{
                "description": "Biochemical Constant; exp(mu_0/RT)/V",
                "units": "per_fmol",
                "symbol": "K",
                "value": 1
            },
            "R":{
                "description":"Universal Gas Constant",
                "units": "J_per_K_mol",
                "symbol": "R",
                "value": 8.31,
            },
            "T":{
                "description": "Temperature",
                "units": "kelvin",
                "symbol": "T",
                "value": 293,
            },
            "q_0":{
                "description":"Molar Quantity",
                "units": "fmol",
                "symbol": "q",
                "value": 1,
            },
          },
          "vars":{
            "e_0": {
                "description": "Generalised Potential",
                "units": "J_per_mol",
                "symbol": "mu",
            },
            "f_0": {
                  "description": "Generalised Flow",
                  "units": "fmol_per_s",
                  "symbol": "v",
            },
          },
          "constitutive_relations":[
            "e_0 - R*T*log(K*q_0)",
          ]
        }

ch_Se_EA = {'componentID': 'ch_Se',
            "vars":{
                "P_0": {
                    "description": "Generalised power",
                    "units": "fJ_per_s",
                    "symbol": "P",
                },
            "P_sum": {
                "description": "Generalised power",
                "units": "fJ_per_s",
                "symbol": "P_sum",
            }
            },
            "state_vars": {
                "A": {
                    "description": "Generalised activity",
                    "units": "fJ",
                    "symbol": "A",
                    "value": 0,
                },
                "E": {
                    "description": "Generalised energy",
                    "units": "fJ",
                    "symbol": "E",
                    "value": 0,
                }
            },
}

m_Se = {
            "description": "Effort Source",
            "domain": "Mechanical",
            "id": "m_Se",
            "metamodel": "SS",
            "ports": {
                "0": {'in':[],'out':[],'type': ' ','direction': 'out'}
            },
            "signals": {
            },
            "params": {
                "e": {
                    "description": "Generalised Potential",
                    "units": "J_per_um",
                    "symbol": "F0",
                    "value": 0
                }
            },
            "vars":{
                "e_0": {
                    "description": "Generalised Potential",
                    "units": "J_per_um",
                    "symbol": "F",
                },
                "f_0": {
                    "description": "Generalised Flow",
                    "units": "um_per_s",
                    "symbol": "dx",
                },
            },
            "constitutive_relations": [
                "e_0 - e"
            ]
        }
m_Se_EA = {'componentID': 'm_Se',
            "vars":{
                "P_0": {
                    "description": "Generalised power",
                    "units": "J_per_s",
                    "symbol": "P",
                },
            "P_sum": {
                "description": "Generalised power",
                "units": "J_per_s",
                "symbol": "P_sum",
            }
            },
            "state_vars": {
                "A": {
                    "description": "Generalised activity",
                    "units": "J",
                    "symbol": "A",
                    "value": 0,
                },
                "E": {
                    "description": "Generalised energy",
                    "units": "J",
                    "symbol": "E",
                    "value": 0,
                }
            },
}
m_Sf = {
            "description": "Flow Source",
            "domain": "Mechanical",
            "id": "m_Sf",
            "metamodel": "SS",
            "ports": {
                "0": {'in':[],'out':[],'type': ' ','direction': 'out'}
            },
            "signals": {
            },
            "params": {
                "f": {
                    "description": "Generalised Flow",
                    "units": "um_per_s",
                    "symbol": "dx0",
                    "value": 0
                }
            },
            "vars":{
                "e_0": {
                    "description": "Generalised Potential",
                    "units": "J_per_um",
                    "symbol": "F",
                },
                "f_0": {
                    "description": "Generalised Flow",
                    "units": "um_per_s",
                    "symbol": "dx",
                },
            },
            "constitutive_relations": [
                "f_0 - f"
            ]
        }
m_Sf_EA = {'componentID': 'm_Sf',
            "vars":{
                "P_0": {
                    "description": "Generalised power",
                    "units": "J_per_s",
                    "symbol": "P",
                },
            "P_sum": {
                "description": "Generalised power",
                "units": "J_per_s",
                "symbol": "P_sum",
            }
            },
            "state_vars": {
                "A": {
                    "description": "Generalised activity",
                    "units": "J",
                    "symbol": "A",
                    "value": 0,
                },
                "E": {
                    "description": "Generalised energy",
                    "units": "J",
                    "symbol": "E",
                    "value": 0,
                }
            },
}
m_R = {
        "description": "Generalised Linear Resistor",
        "domain": "Mechanical",
        "id": "m_R",
        "metamodel": "R",
        "ports": {
            "0": {'in':[],'out':[],'type': ' ','direction': 'in'}
        },
        "signals": {
        },
        "params": {
            "r": {
                "description": "Resistance",
                "value": 1,
                "units": "J_s_per_um2",
                "symbol": "eta",
                "range": [
                    0,
                    "inf"
                ]
            }
        },
        "vars":{
            "e_0": {
                "description": "Generalised Potential",
                "units": "J_per_um",
                "symbol": "F",
            },
            "f_0": {
                "description": "Generalised Flow",
                "units": "um_per_s",
                "symbol": "dx",
            },
        },
        "constitutive_relations": [
            "e_0 - f_0*r"
        ]
    }

m_R_EA = {'componentID': 'm_R',
        "vars":{
            "P_0": {
                "description": "Generalised power",
                "units": "J_per_s",
                "symbol": "P",
            },
            "P_sum": {
                "description": "Generalised power",
                "units": "J_per_s",
                "symbol": "P_sum",
            }
        },
        "state_vars": {
            "A": {
                "description": "Generalised activity",
                "units": "J",
                "symbol": "A",
                "value": 0,
            },
            "E": {
                "description": "Generalised energy",
                "units": "J",
                "symbol": "E",
                "value": 0,
            }
        },
}

m_C = {
        "description": "Generalised Linear Capacitor",
        "domain": "Mechanical",
        "id": "m_C",
        "metamodel": "C",
        "ports": {
            "0": {'in':[],'out':[],'type': ' ','direction': 'in'}
        },
        "signals": {
        },
        "params": {
            "C": {
                "description": "Capacitance",
                "value": 1,
                "units": "um2_per_J",
                "symbol": "mC",
                "range": [
                    0,
                    "inf"
                ]
            },
            "q_init": {
                "description": "Initial generalised position",
                "units": "um",
                "symbol": "x_init",
                "value": 1,
            },
        },
        "vars":{
            "e_0": {
                "description": "Generalised Potential",
                "units": "J_per_um",
                "symbol": "F",
            },
            "f_0": {
                "description": "Generalised Flow",
                "units": "um_per_s",
                "symbol": "dx",
            },
        },
        "state_vars": {
            "q_0": {
                "description": "Generalised position",
                "units": "um",
                "symbol": "x",
                "value": "x_init",
            },
        },
        "constitutive_relations": [
            "q_0 - C * e_0",
            "ode(q_0,t) - f_0"
        ]
    }

m_C_EA = {'componentID': 'm_C',
        "vars":{
            "P_0": {
                "description": "Generalised power",
                "units": "J_per_s",
                "symbol": "P",
            },
            "P_sum": {
                "description": "Generalised power",
                "units": "J_per_s",
                "symbol": "P_sum",
            }
        },
        "state_vars": {
            "A": {
                "description": "Generalised activity",
                "units": "J",
                "symbol": "A",
                "value": 0,
            },
            "E": {
                "description": "Generalised energy",
                "units": "J",
                "symbol": "E",
                "value": 0,
            }
        },
}

cv_R = {
        "description": "Generalised Linear Resistor",
        "domain": "Cardiovascular",
        "id": "cv_R",
        "metamodel": "R",
        "ports": {
            "0": {'in':[],'out':[],'type': ' ','direction': 'in'}
        },
        "signals": {
        },
        "params": {
            "r": {
                "description": "Resistance",
                "value": 1,
                "units": "mmHg_s_per_mL",
                "symbol": "R",
                "range": [
                    0,
                    "inf"
                ]
            }
        },
        "vars":{
            "e_0": {
                "description": "Generalised Potential",
                "units": "mmHg",
                "symbol": "mu",
            },
            "f_0": {
                "description": "Generalised Flow",
                "units": "mL_per_s",
                "symbol": "v",
            },
        },
        "constitutive_relations": [
            "e_0 - f_0*r"
        ]
    }

cv_R_EA = {'componentID': 'cv_R',
        "vars":{
            "P_0": {
                "description": "Generalised power",
                "units": "mmHg_mL_per_s",
                "symbol": "P",
            },
            "P_sum": {
                "description": "Generalised power",
                "units": "mmHg_mL_per_s",
                "symbol": "P_sum",
            }
        },
        "state_vars": {
            "A": {
                "description": "Generalised activity",
                "units": "mmHg_mL",
                "symbol": "A",
                "value": 0,
            },
            "E": {
                "description": "Generalised energy",
                "units": "mmHg_mL",
                "symbol": "E",
                "value": 0,
            }
        },
}

cv_C = {
        "description": "Generalised Linear Capacitor",
        "domain": "Cardiovascular",
        "id": "cv_C",
        "metamodel": "C",
        "ports": {
            "0": {'in':[],'out':[],'type': ' ','direction': 'in'}
        },
        "signals": {
        },
        "params": {
            "C": {
                "description": "Capacitance",
                "value": 1,
                "units": "mL_per_mmHg",
                "symbol": "C",
                "range": [
                    0,
                    "inf"
                ]
            },
            "q_init": {
                "description": "Initial generalised position",
                "units": "mL",
                "symbol": "q_init",
                "value": 1,
            },
        },
        "vars":{
            "e_0": {
                "description": "Generalised Potential",
                "units": "mmHg",
                "symbol": "mu",
            },
            "f_0": {
                "description": "Generalised Flow",
                "units": "mL_per_s",
                "symbol": "v",
            },
        },
        "state_vars": {
            "q_0": {
                "description": "Generalised position",
                "units": "mL",
                "symbol": "q",
                "value": "q_init",
            },
        },
        "constitutive_relations": [
            "q_0 - C * e_0",
            "ode(q_0,t) - f_0"
        ]
    }

cv_C_EA = {'componentID': 'cv_C',
        "vars":{
            "P_0": {
                "description": "Generalised power",
                "units": "mmHg_mL_per_s",
                "symbol": "P",
            },
            "P_sum": {
                "description": "Generalised power",
                "units": "mmHg_mL_per_s",
                "symbol": "P_sum",
            }
        },
        "state_vars": {
            "A": {
                "description": "Generalised activity",
                "units": "mmHg_mL",
                "symbol": "A",
                "value": 0,
            },
            "E": {
                "description": "Generalised energy",
                "units": "mmHg_mL",
                "symbol": "E",
                "value": 0,
            }
        },
}

cv_I = {
        "description": "Generalised Linear Inductor",
        "domain": "Cardiovascular",
        "id": "cv_I",
        "metamodel": "I",
        "ports": {
            "0": {'in':[],'out':[],'type': ' ','direction': 'in'}
        },
        "signals": {
        },
        "vars":{
            "e_0": {
                "description": "Generalised Potential",
                "units": "mmHg",
                "symbol": "mu",
            },
            "f_0": {
                "description": "Generalised Flow",
                "units": "mL_per_s",
                "symbol": "v",
            },
            "P_0": {
                "description": "Generalised power",
                "units": "mmHg_mL_per_s",
                "symbol": "P",
            }
        },
        "state_vars": {
            "p_0": {
                "description": "Generalised momentum",
                "units": "mmHg_mL2_per_s3",
                "symbol": "p",
                "value": "p_init"
            },
             "A": {
                "description": "Generalised activity",
                "units": "mmHg_mL",
                "symbol": "A",
                "value": 0,
            },
            "E": {
                "description": "Generalised energy",
                "units": "mmHg_mL",
                "symbol": "E",
                "value": 0,
            }
        },
        "params": {
            "L": {
                "description": "Inductance",
                "value": 1,
                "units": "mmHg_s2_per_ml",
                "symbol": "L",
                "range": [
                    0,
                    "inf"
                ]
            },
            "p_init": {
                "description": "Initial generalised momentum",
                "units": "mmHg_mL2_per_s3",
                "symbol": "p_init",
                "value": 1
            }   
        },
        "constitutive_relations": [
            "p_0 - L*f_0",
            "ode(p_0, t) - e_0"
        ]
    }

cv_I_EA = {'componentID': 'cv_I',
        "vars":{
            "P_0": {
                "description": "Generalised power",
                "units": "mmHg_mL_per_s",
                "symbol": "P",
            },
            "P_sum": {
                "description": "Generalised power",
                "units": "mmHg_mL_per_s",
                "symbol": "P_sum",
            }
        },
        "state_vars": {
            "A": {
                "description": "Generalised activity",
                "units": "mmHg_mL",
                "symbol": "A",
                "value": 0,
            },
            "E": {
                "description": "Generalised energy",
                "units": "mmHg_mL",
                "symbol": "E",
                "value": 0,
            }
        },
}

def BG_components():
    components = {
        "R": R,
        "C": C,
        "I": I,
        "Se": Se,
        "Sf": Sf,
        "TF": TF,
        "GY": GY,
        "Ce": Ce,
        "Re": Re,
        "zF": zF,
        "ch_Se": ch_Se,
        "m_Se": m_Se,
        "m_Sf": m_Sf,
        "m_R": m_R,
        "m_C": m_C,
        "cv_R": cv_R,
        "cv_C": cv_C,
        "cv_I": cv_I,
    }
    return components 

def BG_EA_components():
    components = {
        "R": R_EA,
        "C": C_EA,
        "I": I_EA,
        "Se": Se_EA,
        "Sf": Sf_EA,
        "TF": TF_EA,
        "Ce": Ce_EA,
        "Re": Re_EA,
        "zF": zF_EA,
        "ch_Se": ch_Se_EA,
        "m_Se": m_Se_EA,
        "m_Sf": m_Sf_EA,
        "m_R": m_R_EA,
        "m_C": m_C_EA,
        "cv_R": cv_R_EA,
        "cv_C": cv_C_EA,
        "cv_I": cv_I_EA,
    }
    return components

if __name__ == "__main__": 

    json_file = 'BG_components.json'
    comp_dict = BG_components()
    with open(json_file, 'w') as f:
        json.dump(comp_dict, f,indent=4)
    json_file = 'BG_EA_components.json'
    comp_dict = BG_EA_components()
    with open(json_file, 'w') as f:
        json.dump(comp_dict, f,indent=4)