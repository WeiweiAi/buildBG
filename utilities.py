import json
def load_json(json_file):
    """
    Load the json file to a dictionary

    Parameters
    ----------
    json_file : str
        The file path of the json file

    Returns
    -------
    comp_dict : dict
        The dictionary of the bond graph model

    """
    with open(json_file) as f:
        comp_dict = json.load(f)
    return comp_dict

def save_json(comp_dict, json_file):
    """
    Save the dictionary to a json file

    Parameters
    ----------
    comp_dict : dict
        The dictionary of the bond graph model
    json_file : str
        The file path of the json file

    Returns
    -------
    None

    side effect
    ------------
    Save the dictionary to a json file

    """
    with open(json_file, 'w') as f:
        json.dump(comp_dict, f,indent=4)