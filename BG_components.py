# Adapted from https://github.com/BondGraphTools/BondGraphTools 
# under the Apache License http://www.apache.org/licenses/LICENSE-2.0
# Add units to the parameters, variables; Add power and activity variables; Add signals; 
# Define individual components
# Note: Natural logarithm ln=log is used in the constitutive relations

# the key of the vars: name_portid, for example, e_0, f_0, s_2, s_3, etc.
# type of ports: defines the type of bond the port, for example, Power, Signal.
# direction: indicating the reference direction of the energy flow, 
# in = energy flows towards the port, out = energy flows away from the port
# '' = direction not defined for signal ports
# causality: indicating the causality of the bond that the port is connected to,
# 'e_in' = effort in, 'e_out' = effort out, '' = causality not defined, 's_in' = signal in, 's_out' = signal out
import json
R = {
        "description": "Generalised Linear Resistor",
        "domain": "Electrical",
        "id": "R",
        "metamodel": "R",
        "ports": {
            "0": {
                "description": "Power port 0",
                "type": "Power",
                'direction': 'in',
                'causality': '',
                'in':[],'out':[]
                },
        },
        "params": {
            "r": {
                "description": "Resistance",
                "value": 1,
                "units": "pOhm",
                "symbol": "R",
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
            "e_0 - f_0*r"
        ]
    }

C = {
        "description": "Generalised Linear Capacitor",
        "domain": "Electrical",
        "id": "C",
        "metamodel": "C",
        "ports": {
            "0": {
                "description": "Power port 0",
                "type": "Power",
                'direction': 'in',
                'causality': '',
                'in':[],'out':[]
                },
        },
        "params": {
            "C": {
                "description": "Capacitance",
                "value": 1,
                "units": "fF",
                "symbol": "C",
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
            "q_0": {
                "description": "Generalised position",
                "units": "fC",
                "symbol": "q",
                "initial value": "q_init"
            },
        },
        "constitutive_relations": [
            "q_0 - C * e_0",
            "ode(q_0,t) - f_0",
        ]
    }

I = {
        "description": "Generalised Linear Inductor",
        "domain": "Electrical",
        "id": "I",
        "metamodel": "I",
        "ports": {
            "0": {
                "description": "Power port 0",
                "type": "Power",
                'direction': 'in',
                'causality': '',
                'in':[],'out':[]
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
            "p_0": {
                "description": "Generalised momentum",
                "units": "volt_s",
                "symbol": "p",
                "initial value": "p_init"
            }
        },
        "params": {
            "L": {
                "description": "Inductance",
                "value": 1,
                "units": "pH",
                "symbol": "L",
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

Se = {
        "description": "Effort Source",
        "domain": "Electrical",
        "id": "Se",
        "metamodel": "SS",
         "ports": {
            "0": {
                "description": "Power port 0",
                "type": "Power",
                'direction': 'out',
                'causality': '',
                'in':[],'out':[]
                },
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


Sf = {
        "description": "Flow Source",
        "domain": "Electrical",
        "id": "Sf",
        "metamodel": "SS",
        "ports": {
            "0": {
                "description": "Power port 0",
                "type": "Power",
                'direction': 'out',
                'causality': '',
                'in':[],'out':[]
                },
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


TF = {
        "description": "Linear Transformer",
        "domain": "Electrical",
        "id": "TF",
        "metamodel": "TF",
        "ports": {
            "0": {
                "description": "Power port 0",
                "type": "Power",
                'direction': 'in',
                'causality': '',
                'in':[],'out':[]
                },
            "1": {
                "description": "Power port 1",
                "type": "Power",
                'direction': 'out',
                'causality': '',
                'in':[],'out':[]
                },
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

GY = {
        "description": "Linear Gyrator",
        "domain": "Electrical-Mechanical",
        "id": "GY",
        "metamodel": "GY",
        "ports": {
            "0": {
                "description": "Power port 0",
                "type": "Power",
                'direction': 'in',
                'causality': '',
                'in':[],'out':[]
                },
            "1": {
                "description": "Power port 1",
                "type": "Power",
                'direction': 'out',
                'causality': '',
                'in':[],'out':[]
                },
        },
        "vars":{
            "e_0": {
                "description": "Generalised Potential",
                "units": "",
                "symbol": "",
            },
            "f_0": {
                "description": "Generalised Flow",
                "units": "",
                "symbol": "",
            },
            "e_1": {
                "description": "Generalised Potential",
                "units": "",
                "symbol": "",
            },
            "f_1": {
                "description": "Generalised Flow",
                "units": "",
                "symbol": "",
            },

        },
        "params": {
            "r": {
                "description": "Ratio",
                "value": 1,
                "symbol": "n",
                "units": "",
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
        "ports": {
            "0": {
                "description": "Power port 0",
                "type": "Power",
                'direction': 'in',
                'causality': '',
                'in':[],'out':[]
                },
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
          },
          "q_0":{
              "description":"Molar Quantity",
              "units": "fmol",
              "symbol": "q",
              "initial value": "q_init",
          },      
        },
        "constitutive_relations":[
          "e_0 - R*T*log(K*q_0)",
          "ode(q_0,t) - f_0"
        ]
     }


Re ={
        "description": "Biochemical Reaction",
        "domain": "biochemical",
        "id": "Re",
        "metamodel":"R",
        "ports": {
            "0": {
                "description": "Power port 0",
                "type": "Power",
                'direction': 'in',
                'causality': '',
                'in':[],'out':[]
                },
            "1": {
                "description": "Power port 1",
                "type": "Power",
                'direction': 'out',
                'causality': '',
                'in':[],'out':[]
                },
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
              },
        },
        "constitutive_relations":[
          "f_0 - kappa*(exp(e_0/R/T) - exp(e_1/R/T))"
        ]
    }

zF = {
        "description": "Biochemical Transformer",
        "domain": "biochemical",
        "id": "zF",
        "metamodel": "TF",
        "ports": {
            "0": {
                "description": "Power port 0",
                "type": "Power",
                'direction': 'in',
                'causality': '',
                'in':[],'out':[]
                },
            "1": {
                "description": "Power port 1",
                "type": "Power",
                'direction': 'out',
                'causality': '',
                'in':[],'out':[]
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

Re_GHK ={
        "description": "Voltage modulated Biochemical Reaction-Goldman-Hodgkin-Katz (GHK) ion channel",
        "domain": "biochemical",
        "id": "Re_GHK",
        "metamodel":"R",
        "ports": {
            "0": {
                "description": "Power port 0",
                "type": "Power",
                'direction': 'in',
                'causality': '',
                'in':[],'out':[]
                },
            "1": {
                "description": "Power port 1",
                "type": "Power",
                'direction': 'out',
                'causality': '',
                'in':[],'out':[]
                },
            "m": {
                "description": "modulation port",
                "type": "Signal",
                'direction': '',
                'causality': 's_in',
                'in':[],'out':[]
                },
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
              },
              "s_m": {
                  "description": "Generalised Potential",
                  "units": "J_per_mol",
                  "symbol": "Am",
              },
        },
        "constitutive_relations":[
          "f_0- Piecewise(( kappa*(exp(e_0/R/T) - exp(e_1/R/T)),s_m==0),(kappa*(s_m/R/T/(exp(s_m/R/T)))*((exp(e_0/R/T) - exp(e_1/R/T))),s_m!=0))"
        ]
    }

ch_Se = {
          "description":"Concentration of Chemical Species",
          "domain": "biochemical",
          "id":"ch_Se",
          "metamodel":"SS",
          "ports": {
            "0": {
                "description": "Power port 0",
                "type": "Power",
                'direction': 'out',
                'causality': '',
                'in':[],'out':[]
                },
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

m_Se = {
            "description": "Effort Source",
            "domain": "Mechanical",
            "id": "m_Se",
            "metamodel": "SS",
            "ports": {
                "0": {
                   "description": "Power port 0",
                   "type": "Power",
                   'direction': 'out',
                   'causality': '',
                   'in':[],'out':[]
                   },
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
m_Sf = {
            "description": "Flow Source",
            "domain": "Mechanical",
            "id": "m_Sf",
            "metamodel": "SS",
            "ports": {
             "0": {
                   "description": "Power port 0",
                   "type": "Power",
                   'direction': 'out',
                   'causality': '',
                   'in':[],'out':[]
                   },
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

m_R = {
        "description": "Generalised Linear Resistor",
        "domain": "Mechanical",
        "id": "m_R",
        "metamodel": "R",
        "ports": {
             "0": {
                   "description": "Power port 0",
                   "type": "Power",
                   'direction': 'in',
                   'causality': '',
                   'in':[],'out':[]
                   },
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

m_C = {
        "description": "Generalised Linear Capacitor",
        "domain": "Mechanical",
        "id": "m_C",
        "metamodel": "C",
        "ports": {
             "0": {
                   "description": "Power port 0",
                   "type": "Power",
                   'direction': 'in',
                   'causality': '',
                   'in':[],'out':[]
                   },
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
            "q_0": {
                "description": "Generalised position",
                "units": "um",
                "symbol": "x",
                "initial value": "x_init",
            },
        },
        "constitutive_relations": [
            "q_0 - C * e_0",
            "ode(q_0,t) - f_0"
        ]
    }


cv_R = {
        "description": "Generalised Linear Resistor",
        "domain": "Cardiovascular",
        "id": "cv_R",
        "metamodel": "R",
        "ports": {
             "0": {
                   "description": "Power port 0",
                   "type": "Power",
                   'direction': 'in',
                   'causality': '',
                   'in':[],'out':[]
                   },
        },
        "params": {
            "r": {
                "description": "Resistance",
                "value": 1,
                "units": "mmHg_s_per_mL",
                "symbol": "R",
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

cv_C = {
        "description": "Generalised Linear Capacitor",
        "domain": "Cardiovascular",
        "id": "cv_C",
        "metamodel": "C",
        "ports": {
             "0": {
                   "description": "Power port 0",
                   "type": "Power",
                   'direction': 'in',
                   'causality': '',
                   'in':[],'out':[]
                   },
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
            "q_0": {
                "description": "Generalised position",
                "units": "mL",
                "symbol": "q",
                "initial value": "q_init",
            },
        },
        "constitutive_relations": [
            "q_0 - C * e_0",
            "ode(q_0,t) - f_0"
        ]
    }

cv_I = {
        "description": "Generalised Linear Inductor",
        "domain": "Cardiovascular",
        "id": "cv_I",
        "metamodel": "I",
        "ports": {
             "0": {
                   "description": "Power port 0",
                   "type": "Power",
                   'direction': 'in',
                   'causality': '',
                   'in':[],'out':[]
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
            "p_0": {
                "description": "Generalised momentum",
                "units": "mmHg_mL2_per_s3",
                "symbol": "p",
                "initial value": "p_init"
            },
        },
        "params": {
            "L": {
                "description": "Inductance",
                "value": 1,
                "units": "mmHg_s2_per_ml",
                "symbol": "L",
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
        "Re_GHK": Re_GHK,
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

if __name__ == "__main__": 

    json_file = 'BG_components.json'
    comp_dict = BG_components()
    with open(json_file, 'w') as f:
        json.dump(comp_dict, f,indent=4)
