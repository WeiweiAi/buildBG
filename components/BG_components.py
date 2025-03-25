# Adapted from https://github.com/BondGraphTools/BondGraphTools 
# under the Apache License http://www.apache.org/licenses/LICENSE-2.0
import json

if __name__ == "__main__": 

    json_file = 'BG_components.json'
    
    with open("BG_biochemical.json") as j1, open("BG_electrical.json") as j2,  open("BG_mechanical.json") as j3, open("BG_hydraulic.json") as j4, open(json_file, 'w') as j5:
        json.dump({1:json.load(j1), 2:json.load(j2), 3:json.load(j3),4:json.load(j4)},j5,indent=4)