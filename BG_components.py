# Adapted from https://github.com/BondGraphTools/BondGraphTools 
# under the Apache License http://www.apache.org/licenses/LICENSE-2.0
# Add units to the parameters, variables;
import copy
def _base_components():
    base_components = {
        "description": "Basic Components",
        "id": "base",
        "components": {
            "0": {
                "description": "Equal Effort Node",
                "metamodel": "0",
                "class": "EqualEffort"
            },
            "1": {
                "description": "Equal Flow Node",
                "class": "EqualFlow",
                "metamodel": "1"
            },
            "R": {
                "description": "Generalised Linear Resistor",
                "class": "Component",
                "metamodel": "R",
                "ports": {
                    "0": {'in':[],'out':[]}
                },
                "params": {
                    "r": {
                        "description": "Resistance",
                        "value": 1,
                        "range": [
                            0,
                            "inf"
                        ]
                    }
                },
                "vars":{
                    "e_0": {
                        "description": "Generalised Potential",
                    },
                    "f_0": {
                        "description": "Generalised Flow",
                    }
                },
                "constitutive_relations": [
                    "e_0 - f_0*r"
                ]
            },
            "C": {
                "description": "Generalised Linear Capacitor",
                "class": "Component",
                "metamodel": "C",
                "ports": {
                    "0": {'in':[],'out':[]}
                },
                "params": {
                    "C": {
                        "description": "Capacitance",
                        "value": 1,
                        "range": [
                            0,
                            "inf"
                        ]
                    },
                    "q_init": {"description": "Generalised position",
                               "value": 1,
                    },
                },
                "vars":{
                    "e_0": {
                        "description": "Generalised Potential",
                    },
                    "f_0": {
                        "description": "Generalised Flow",
                    }
                },
                "state_vars": {
                    "q_0": {
                        "description": "Generalised position",
                        "value": "q_init",
                    }
                },
                "constitutive_relations": [
                    "q_0 - C * e_0",
                    "dq_0 - f_0"
                ]
            },
            "I": {
                "description": "Generalised Linear Inductor",
                "class": "Component",
                "metamodel": "I",
                "ports": {
                    "0": {'in':[],'out':[]}
                },
                "vars":{
                    "e_0": {
                        "description": "Generalised Potential",
                    },
                    "f_0": {
                        "description": "Generalised Flow",
                    }
                },
                "state_vars": {
                    "p_0": {
                        "description": "Generalised momentum"
                    },
                },
                "params": {
                    "L": {
                        "description": "Inductance",
                        "value": 1,
                        "range": [
                            0,
                            "inf"
                        ]
                    }
                },
                "constitutive_relations": [
                    "p_0 - L*f_0",
                    "dp_0 - e_0"
                ]
            },
            "Se": {
                "description": "Effort Source",
                "class": "Component",
                "metamodel": "SS",
                "ports": {
                    "0": {'in':[],'out':[]}
                },
                "params": {
                    "e": {
                        "description": "Generalised Potential",
                        "value": 0
                    }
                },
                "vars":{
                    "e_0": {
                        "description": "Generalised Potential",
                    }
                },
                "constitutive_relations": [
                    "e_0 = e"
                ]
            },
            "Sf": {
                "description": "Flow Source",
                "class": "Component",
                "metamodel": "SS",
                "ports": {
                    "0": {'in':[],'out':[]}
                },
                "params": {
                    "f": {
                        "description": "Generalised Flow",
                        "value": 0
                    }
                },
                "vars":{
                    "f_0": {
                        "description": "Generalised Flow",
                    }
                },
                "constitutive_relations": [
                    "f_0 = f"
                ]
            },
            "TF": {
                "description": "Linear Transformer",
                "class": "Component",
                "metamodel": "TF",
                "ports": {
                    "0": {
                        "description": "Primary",
                        'in':[],'out':[],
                    },
                    "1": {
                        "description": "Secondary",
                        'in':[],'out':[]
                    }
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
                    "e_1 - r * e_0",
                    "f_0 + r * f_1"
                ]
            },
            "GY": {
                "description": "Linear Gyrator",
                "class": "Component",
                "metamodel": "GY",
                "ports": {
                    "0": {
                        "description": "Primary",
                        'in':[],'out':[]
                    },
                    "1": {
                        "description": "Secondary",
                        'in':[],'out':[]
                    }
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
                    "e_1 + r*f_0",
                    "e_0 - r*f_1"
                ]
            },
            "SS": {
                "description": "Source Sensor",
                "class": "Component",
                "metamodel": "SS",
                "ports": {
                    "0": {'in':[],'out':[]}
                },
                "params": {
                    "f": {},
                    "e": {}
                },
                "constitutive_relations": ["e_0 - e", "f_0+f"]
            },
            "PH": {
                "class": "PortHamiltonian",
                "metamodel": "PH",
                "description": "Port Hamiltonian",
                "params": {}
            }
        }
    }
    return base_components

def _biochem_components():

    biochem_components= {
      "description": "Biochemical Components",
      "id":"BioChem",
      "components":{
        "Ce":{
          "description":"Concentration of Chemical Species",
          "class": "Component",
          "metamodel":"C",
          "ports":{
            "0":{'in':[],'out':[]},
          },
          "params":{
            "K":{"description": "Biochemical Constant; exp(mu_0/RT)/V",
                "value": 1,},
            "R":{"description":"Universal Gas Constant",
                 "value": 8.31,},
            "T":{"description": "Temperature",
                 "value": 293,},
            "q_init":{"description": "Initial Molar Quantity",
                    "value": 1,}
          },
          "state_vars":{ 
            "q_0":{"description":"Molar Quantity",
                   "value": "q_init",},

          },
          "vars":{
            "e_0": {
                "description": "Generalised Potential",
            },
            "f_0": {
                "description": "Generalised Flow ",
            },
          },
          "constitutive_relations":[
            "e_0 = R*T*log(K*q_0)",
            "ode(q_0,t) = f_0"
          ]
        },
        "Re":{
          "description": "Biochemical Reaction",
          "class": "SymmetricComponent",
          "metamodel":"R",
          "ports":{
            "0":{'in':[],'out':[]},
            "1":{'in':[],'out':[]}
          },
          "params":{
            "kappa":{"description":"Reaction Rate",
                     "value": 1,},
            "R":{"description":"Universal Gas Constant",
                 "value": 8.31,},
            "T":{"description": "Temperature",
                 "value": 293,}
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
                    }
                },
          "constitutive_relations":[
            "f_0 = kappa*(exp(e_0/R/T) - exp(e_1/R/T))"
          ]
        },
      "Y":{
        "description": "Stoichiometry",
        "class": "MixedPort",
        "metamodel": "TF",
        "ports":{
          "0": {
            "description":"Affinity",
            "r": -1,
            'in':[],'out':[]
          },
          "i": {
            "description":"Chemical Power",
            "r": 1,
            'in':[],'out':[]
          }
        },
        "constitutive_relations":[
          "sum(r_i*e_i)",
          "f_0/r_0 - f_i/r_i"
        ]
      }
      }
    }
    return biochem_components

def e_components_units():
    e_components = copy.deepcopy(_base_components())
    e_components['components']['e_Se']= copy.deepcopy(_base_components()['components']['Se'])
    e_components['components']['e_Sf']= copy.deepcopy(_base_components()['components']['Sf'])
    # remove the original Se, Sf from the base_components
    e_components['components'].pop('Se')
    e_components['components'].pop('Sf')
    e_components['components']['e_Se']['params']['e']['units'] = 'volt'
    e_components['components']['e_Sf']['params']['f']['units'] = 'fA'
    e_components['components']['e_Se']['params']['e']['symbol'] = 'V0'
    e_components['components']['e_Sf']['params']['f']['symbol'] = 'I0'
    e_components['components']['e_Sf']['vars']['f_0']['units'] = 'fA'
    e_components['components']['e_Se']['vars']['e_0']['units'] = 'volt'
    e_components['components']['e_Sf']['vars']['f_0']['symbol'] = 'I'
    e_components['components']['e_Se']['vars']['e_0']['symbol'] = 'V'
    e_components['components']['e_Sf']['vars']['f_0']['IOType'] = 'out'
    e_components['components']['e_Se']['vars']['e_0']['IOType'] = 'out'

    e_components['components']['C']['params']['C']['units'] = 'fF'
    e_components['components']['C']['params']['q_init']['units'] = 'fC'
    e_components['components']['C']['params']['C']['symbol'] = 'C'
    e_components['components']['C']['params']['q_init']['symbol'] = 'q_init'

    e_components['components']['C']['vars']['e_0']['units'] = 'volt'
    e_components['components']['C']['vars']['f_0']['units'] = 'fA'
    e_components['components']['C']['vars']['e_0']['symbol'] = 'E'
    e_components['components']['C']['vars']['f_0']['symbol'] = 'I'
    e_components['components']['C']['vars']['e_0']['IOType'] = 'out'
    e_components['components']['C']['vars']['f_0']['IOType'] = 'in'
    e_components['components']['C']['state_vars']['q_0']['units'] = 'fC'
    e_components['components']['C']['state_vars']['q_0']['symbol'] = 'q'
    e_components['components']['C']['state_vars']['q_0']['IOType'] = 'out'

    e_components['components']['R']['params']['r']['units'] = 'pOhm'
    e_components['components']['R']['params']['r']['symbol'] = 'R'

    e_components['components']['R']['vars']['e_0']['units'] = 'volt'
    e_components['components']['R']['vars']['f_0']['units'] = 'fA'
    e_components['components']['R']['vars']['e_0']['symbol'] = 'E'
    e_components['components']['R']['vars']['f_0']['symbol'] = 'I'
    e_components['components']['R']['vars']['e_0']['IOType'] = 'in'
    e_components['components']['R']['vars']['f_0']['IOType'] = 'out'
    
   
    e_components['components']['zF'] = copy.deepcopy(_base_components()['components']['TF'])
    e_components['components']['zF']['params']['r']['units'] = 'dimensionless'
    e_components['components']['zF']['params']['r']['symbol'] = 'z'
    # add F to the parameters in the zF components for physiology applications
    e_components['components']['zF']['params']['F'] = {}
    e_components['components']['zF']['params']['F']['units'] = 'C_per_mol'
    e_components['components']['zF']['params']['F']['symbol'] = 'F'

    e_components['components']['zF']['vars']['e_0']['units'] = 'volt'
    e_components['components']['zF']['vars']['f_0']['units'] = 'fA' 
    e_components['components']['zF']['vars']['e_1']['units'] = 'J_per_mol'
    e_components['components']['zF']['vars']['f_1']['units'] = 'fmol_per_s' 
    e_components['components']['zF']['vars']['e_0']['symbol'] = 'V'
    e_components['components']['zF']['vars']['f_0']['symbol'] = 'I' 
    e_components['components']['zF']['vars']['e_1']['symbol'] = 'mu'
    e_components['components']['zF']['vars']['f_1']['symbol'] = 'v'    
    e_components['components']['zF']['vars']['e_0']['IOType'] = 'in'
    e_components['components']['zF']['vars']['f_0']['IOType'] = 'out' 
    e_components['components']['zF']['vars']['e_1']['IOType'] = 'out'
    e_components['components']['zF']['vars']['f_1']['IOType'] = 'in' 

    return e_components

def biochem_components_units():
    # add components Se, Sf, TF from base_components
    biochem_components = copy.deepcopy(_biochem_components())
    biochem_components['components']['Se'] = copy.deepcopy(_base_components()['components']['Se'])
    biochem_components['components']['Sf'] = copy.deepcopy(_base_components()['components']['Sf'])
    biochem_components['components']['TF'] = copy.deepcopy(_base_components()['components']['TF'])
    # Add units to the parameters in the biochem components for physiology applications    
    biochem_components['components']['Ce']['params']['K']['units'] = 'per_fmol'
    biochem_components['components']['Ce']['params']['R']['units'] = 'J_per_K_mol'
    biochem_components['components']['Ce']['params']['T']['units'] = 'kelvin'
    biochem_components['components']['Ce']['params']['q_init']['units'] = 'fmol'    
    biochem_components['components']['Ce']['params']['K']['symbol'] = 'K'
    biochem_components['components']['Ce']['params']['R']['symbol'] = 'R'
    biochem_components['components']['Ce']['params']['T']['symbol'] = 'T'
    biochem_components['components']['Ce']['params']['q_init']['symbol'] = 'q_init'

    biochem_components['components']['Ce']['vars']['e_0']['units'] = 'J_per_mol'
    biochem_components['components']['Ce']['vars']['f_0']['units'] = 'fmol_per_s'
    biochem_components['components']['Ce']['vars']['e_0']['symbol'] = 'mu'
    biochem_components['components']['Ce']['vars']['f_0']['symbol'] = 'v'
    biochem_components['components']['Ce']['vars']['e_0']['IOType'] = 'out'
    biochem_components['components']['Ce']['vars']['f_0']['IOType'] = 'in'
    biochem_components['components']['Ce']['state_vars']['q_0']['units'] = 'fmol'
    biochem_components['components']['Ce']['state_vars']['q_0']['symbol'] = 'q'
    biochem_components['components']['Ce']['state_vars']['q_0']['IOType'] = 'out'

    biochem_components['components']['ch_Se'] = copy.deepcopy(biochem_components['components']['Ce'])
    biochem_components['components']['ch_Se']['params']['q_init']['symbol'] = 'q'
    biochem_components['components']['ch_Se'].pop('state_vars')
    biochem_components['components']['ch_Se']['vars'].pop('f_0')

    biochem_components['components']['Re']['params']['kappa']['units'] = 'fmol_per_s'
    biochem_components['components']['Re']['params']['R']['units'] = 'J_per_K_mol'
    biochem_components['components']['Re']['params']['T']['units'] = 'kelvin'
    biochem_components['components']['Re']['params']['kappa']['symbol'] = 'kappa'
    biochem_components['components']['Re']['params']['R']['symbol'] = 'R'
    biochem_components['components']['Re']['params']['T']['symbol'] = 'T'

    biochem_components['components']['Re']['vars']['e_0']['units'] = 'J_per_mol'
    biochem_components['components']['Re']['vars']['f_0']['units'] = 'fmol_per_s'
    biochem_components['components']['Re']['vars']['e_1']['units'] = 'J_per_mol'
    biochem_components['components']['Re']['vars']['e_0']['symbol'] = 'A_f'
    biochem_components['components']['Re']['vars']['f_0']['symbol'] = 'v'
    biochem_components['components']['Re']['vars']['e_1']['symbol'] = 'A_r'
    biochem_components['components']['Re']['vars']['e_0']['IOType'] = 'in'
    biochem_components['components']['Re']['vars']['f_0']['IOType'] = 'out'
    biochem_components['components']['Re']['vars']['e_1']['IOType'] = 'in'

    
    biochem_components['components']['Se']['params']['e']['units'] = 'J_per_mol'
    biochem_components['components']['Sf']['params']['f']['units'] = 'fmol_per_s'
    biochem_components['components']['Se']['params']['e']['symbol'] = 'mu0'
    biochem_components['components']['Sf']['params']['f']['symbol'] = 'v0'
    
    biochem_components['components']['Se']['vars']['e_0']['units'] = 'J_per_mol'
    biochem_components['components']['Sf']['vars']['f_0']['units'] = 'fmol_per_s'
    biochem_components['components']['Se']['vars']['e_0']['symbol'] = 'mu'
    biochem_components['components']['Sf']['vars']['f_0']['symbol'] = 'v'
    biochem_components['components']['Se']['vars']['e_0']['IOType'] = 'out'
    biochem_components['components']['Sf']['vars']['f_0']['IOType'] = 'out'


    biochem_components['components']['TF']['params']['r']['units'] = 'dimensionless'
    biochem_components['components']['TF']['params']['r']['symbol'] = 'n' 
   
    
    biochem_components['components']['TF']['vars']['e_0']['units'] = 'J_per_mol'
    biochem_components['components']['TF']['vars']['f_0']['units'] = 'fmol_per_s'
    biochem_components['components']['TF']['vars']['e_1']['units'] = 'J_per_mol'
    biochem_components['components']['TF']['vars']['f_1']['units'] = 'fmol_per_s'      
    biochem_components['components']['TF']['vars']['e_0']['symbol'] = 'mu_in'
    biochem_components['components']['TF']['vars']['f_0']['symbol'] = 'v_out'
    biochem_components['components']['TF']['vars']['e_1']['symbol'] = 'mu_out'
    biochem_components['components']['TF']['vars']['f_1']['symbol'] = 'v_in'
    biochem_components['components']['TF']['vars']['e_0']['IOType'] = 'in'
    biochem_components['components']['TF']['vars']['f_0']['IOType'] = 'out'
    biochem_components['components']['TF']['vars']['e_1']['IOType'] = 'out'
    biochem_components['components']['TF']['vars']['f_1']['IOType'] = 'in'
    

    return biochem_components

def m_components_units(): # for mechanics
    m_components = copy.deepcopy(_base_components())
    m_components['components']['m_Se']= copy.deepcopy(_base_components()['components']['Se'])
    m_components['components']['m_Sf']= copy.deepcopy(_base_components()['components']['Sf'])
    # remove the original Se, Sf from the base_components
    m_components['components'].pop('Se')
    m_components['components'].pop('Sf')
    m_components['components']['m_Se']['params']['e']['units'] = 'J_per_um'
    m_components['components']['m_Sf']['params']['f']['units'] = 'um_per_s'
    m_components['components']['m_Se']['params']['e']['symbol'] = 'F0'
    m_components['components']['m_Sf']['params']['f']['symbol'] = 'dx0'
    m_components['components']['m_Sf']['vars']['f_0']['units'] = 'um_per_s'
    m_components['components']['m_Se']['vars']['e_0']['units'] = 'J_per_um'
    m_components['components']['m_Sf']['vars']['f_0']['symbol'] = 'dx'
    m_components['components']['m_Se']['vars']['e_0']['symbol'] = 'F'
    m_components['components']['m_Sf']['vars']['f_0']['IOType'] = 'out'
    m_components['components']['m_Se']['vars']['e_0']['IOType'] = 'out'

    m_components['components']['m_R']= copy.deepcopy(_base_components()['components']['R'])
    m_components['components']['m_C']= copy.deepcopy(_base_components()['components']['C'])
    # remove the original Se, Sf from the base_components
    m_components['components'].pop('R')
    m_components['components'].pop('C')

    m_components['components']['m_R']['params']['r']['units'] = 'J_s_per_um2'
    m_components['components']['m_R']['params']['r']['symbol'] = 'eta'

    m_components['components']['m_R']['vars']['e_0']['units'] = 'J_per_um'
    m_components['components']['m_R']['vars']['f_0']['units'] = 'um_per_s'
    m_components['components']['m_R']['vars']['e_0']['symbol'] = 'F'
    m_components['components']['m_R']['vars']['f_0']['symbol'] = 'dx'
    m_components['components']['m_R']['vars']['e_0']['IOType'] = 'in'
    m_components['components']['m_R']['vars']['f_0']['IOType'] = 'out'

    m_components['components']['m_C']['params']['C']['units'] = 'um2_per_J'
    m_components['components']['m_C']['params']['q_init']['units'] = 'um'
    m_components['components']['m_C']['params']['C']['symbol'] = 'mC'
    m_components['components']['m_C']['params']['q_init']['symbol'] = 'x_init'

    m_components['components']['m_C']['vars']['e_0']['units'] = 'J_per_um'
    m_components['components']['m_C']['vars']['f_0']['units'] = 'um_per_s'
    m_components['components']['m_C']['vars']['e_0']['symbol'] = 'F'
    m_components['components']['m_C']['vars']['f_0']['symbol'] = 'dx'
    m_components['components']['m_C']['vars']['e_0']['IOType'] = 'out'
    m_components['components']['m_C']['vars']['f_0']['IOType'] = 'in'
    m_components['components']['m_C']['state_vars']['q_0']['units'] = 'um'
    m_components['components']['m_C']['state_vars']['q_0']['symbol'] = 'x'
    m_components['components']['m_C']['state_vars']['q_0']['IOType'] = 'out'
   
    

    return m_components    