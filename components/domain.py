# Adapted from https://github.com/BondGraphTools/BondGraphTools 
# under the Apache License http://www.apache.org/licenses/LICENSE-2.0
# Add units to the parameters, variables; 
# Adopted Table 1.2 from  Borutzky, Wolfgang. Bond graph modelling of engineering systems. Vol. 103. New York: springer, 2011.
# Note: Natural logarithm ln=log is used in the constitutive relations

import json
BondElement_domain = {'e':'electrical', 'm':'mechanical', 'hyd':'hydraulic', 'bio':'biochemical','bio-e':'biochemical-electrical'}

if __name__ == "__main__": 

    json_file = 'BG_domain.json'
    with open(json_file, 'w') as f:
        json.dump(BondElement_domain, f,indent=4)
