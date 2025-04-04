import json
from sympy import *
import pandas as pd
from scipy import integrate
import numpy as np
from pathlib import Path
import csv
import libsbml
from networkx.readwrite import json_graph
from pyvis.network import Network
import webbrowser
import os

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
        The dictionary of the json file

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

def load_matrix(matrix):
    """
    Load stoichiometric matrix from a csv file
    The components are defined in BG_components.
    The ID of the components is the key of the dictionary BG_components.
    The components in the row take the potential as the input and the flow as the output, and we call them F components.
    The components in the column take the flow as the input and the potential as the output, and we call them E components.
    The csv file should have the following format (* means blank):
    *       *     *      fName fName
    *       *     *       fID   fID
    *       *     *      fPort fPort
    eName  eID  ePort      0    1 
    eName  eID  ePort      1    0
    
    Parameters
    ----------
    matrix : str
        The file path of the stoichiometric matrix
    
    Returns
    -------
    eName : list
        A list of E component (e_out, f_in) names
    eID : list
        A list of E component (e_out, f_in) IDs
    ePort : list
        A list of E component (e_out, f_in) port number
    fName : list
        A list of F component (e_in, f_out) names
    fID : list
        A list of F component (e_in, f_out) IDs
    fPort : list
        A list of F component (e_in, f_out) port number
    N : numpy.ndarray
        The stoichiometric matrix
    """
    startC=3
    N = []
    eName=[]
    eID=[]
    ePort=[]
    with open(matrix,'r') as f:
        reader = csv.reader(f,delimiter=',')
        line_count = 0
        for row in reader:
            if line_count ==0:
                fName=row[startC:]
                line_count += 1
            elif line_count ==1:
                fID=row[startC:]
                line_count += 1
            elif line_count == 2:
                fPort=row[startC:]
                line_count += 1
            else:
                N.append(row[startC:])
                eName.append(row[0])
                eID.append(row[1])
                ePort.append(row[2])
        f.close()
    
    return eName, eID, ePort, fName, fID, fPort, np.array(N).astype(int)

def infix_to_mathml(ode_var,infix, voi='',version='1.1'):
    """
    Convert the infix string to mathML string defined in CellML specification

    Parameters
    ----------
    ode_var : str
        The derivative variable name
    infix : str
        The infix string to be converted to mathML string
    voi : str, optional
        The variable of integration. The default is ''.
    version : str, optional
        The version of the CellML specification. The default is '1.1'.

    Returns
    -------
    str
        The mathML string defined in CellML specification

    """

    if voi!='':
        preforumla = '<apply> <eq/> <apply> <diff/> <bvar> <ci>'+ voi + '</ci> </bvar> <ci>' + ode_var + '</ci> </apply> '
    else:
        preforumla = '<apply> <eq/> <ci>'+ ode_var + '</ci>'    
    postformula = ' </apply> '
    # replace log to ln in infix string
    infix = infix.replace('log', 'ln')
    # replace fabs to abs in infix string
    infix = infix.replace('fabs', 'abs')
    p = libsbml.parseL3Formula (infix)
    mathstr = libsbml.writeMathMLToString (p)
     # remove the <math> tags in the mathML string, and the namespace declaration will be added later according to the CellML specification
    mathstr = mathstr.replace ('<math xmlns="http://www.w3.org/1998/Math/MathML">', '')
    mathstr = mathstr.replace ('</math>', '')
    mathstr = mathstr.replace ('<?xml version="1.0" encoding="UTF-8"?>', '')
    # temporary solution to add cellml units for constant in the mathML string, replace <cn type="integer"> to <cn cellml:units="dimensionless">
    # check if <cn type="integer"> is in the mathML string
    if '<cn type="integer">' in mathstr or '<cn type="real">' in mathstr:
        mathstr = mathstr.replace ('<cn type="integer">', '<cn cellml:units="dimensionless">')
        mathstr = mathstr.replace ('<cn type="real">', '<cn cellml:units="dimensionless">')
        # add left side of the equation       
    mathstr = preforumla + mathstr + postformula
    # add the cellml namespace to the mathML string
    if '<cn cellml:units="dimensionless">' in mathstr:
        mathstr = f'<apply xmlns:cellml="http://www.cellml.org/cellml/{version}#">' + mathstr + ' </apply>'
    else:
        mathstr = '<apply>' + mathstr + ' </apply>'
    return mathstr

def save_nxBG_json(G, filename='nx_BG.json'):
    """
    Save the bond graph in json format using networkx and json_graph
    Parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx
    filename : str, optional
        The name of the file to save the bond graph. The default is 'nx_BG.json'.
    """
    with open(filename, 'w') as f:
        json.dump(json_graph.node_link_data(G,edges="edges"), f, indent=4)

def load_nxBG_json(filename):
    """
    Load the bond graph in json format using networkx and json_graph
    Parameters
    ----------
    filename : str
        The name of the file to load the bond graph.
    Returns
    -------
    G : nx.DiGraph
        The bond graph built using networkx
    """

    with open(filename, 'r') as f:
        data = json.load(f)
    G = json_graph.node_link_graph(data, edges="edges", directed=True)
    return G

def _pyvis_BG(G):
    """
    Convert the bond graph to pyvis network for visualization

    Parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx

    Returns
    -------
    pyvis_net : pyvis.network.Network
        The bond graph in pyvis network format
    """
    pyvis_net = Network(cdn_resources='in_line', directed=True)
    pyvis_net.from_nx(G)
    return pyvis_net

def _pyvis_net_fomat(pyvis_net):
    """
    Format the bond graph (pyvis_net) for visualization
        
    """ 
    for node in pyvis_net.nodes:
        if 'a' in node.keys():
            if node['a']=='BondElement':
                node['shape']='box'
                node['size']=20
                node['color']={'background':'cyan'}
                node['group']='BondElement'
            elif node['a']=='JunctionStructure':
                node['size']=5
                if node['subClass']=='ZeroJunctionStructure':
                    node['shape']='dot'
                    node['color']={'background':'blue'}
                elif node['subClass']=='OneJunctionStructure':
                    node['shape']='dot'
                    node['color']={'background':'black'}
                elif node['subClass']=='TF' or node['subClass']=='MTF':
                    node['shape']='square'
                    node['color']={'background':'yellow'}
                elif node['subClass']=='GY' or node['subClass']=='MGY':
                    node['shape']='square'
                    node['color']={'background':'blue'}
            elif node['a']=='PowerPort':
                node['shape']='diamond'
                node['size']=5
                node['color']={'background':'white'}
                node['group']='PowerPort'
            elif node['a']=='SignalPort':
                node['shape']='hexagon'
                node['size']=5
                node['color']={'background':'black'}
            node['label']=node['id']
            node['font']={'size':10} 
    for edge in pyvis_net.edges:
        if 'a' in edge.keys():
            if edge['a']=='PowerBond':
                edge['color']={'color':'blue'}
                edge['length']=100
                edge['width']=2
            elif edge['a']=='SignalBond':
                edge['color']={'color':'black'}
        else:
            edge['dashes']='true'
            edge['color']={'color':'gray'}
            edge['length']=10
            edge['width']=1

def _nxBG_pyvis_net(pyvis_net, filename='nx_BG.html'):
    """
    Save the bond graph in pyvis network format to html for visualization

    Parameters
    ----------
    pyvis_net : pyvis.network.Network
        The bond graph in pyvis network format
    filename : str, optional
        The name of the file to save the bond graph. The default is 'nx_BG.html'.

    Returns
    -------
    None
        
    """
    pyvis_net.set_options("""
        var options = {
            "edges": {
                "smooth": true,
                "width": 2
            },
            "physics": {
                "enabled": true,
                "barnesHut": {
                  "gravitationalConstant": -20000,
                  "centralGravity": 0.8,
                  "springLength": 100,  
                  "springConstant": 0.8
                },
                "repulsion": {
                  "nodeDistance": 20, 
                  "springLength": 20    
                },
                "hierarchicalRepulsion": {
                  "nodeDistance": 10,
                  "springLength": 5 
                },
                "minVelocity": 0.75
            }
        }
    """)
   # Generate the HTML content
    html_content = pyvis_net.generate_html()

    # Remove the div with class "card"
    html_content = html_content.replace('<div class="card" style="width: 100%">', '').replace('<div id="mynetwork" class="card-body">', '<div id="mynetwork">').replace('</div>', '', 1)

    # Ensure the graph takes the whole page
    html_content = html_content.replace('<div id="mynetwork">', '<div id="mynetwork" style="width: 100%; height: 100vh;">')
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(html_content)

def save_nxBG_html(G, filename='nx_BG.html'):
    """
    Save the bond graph in pyvis network format to html for visualization

    Parameters
    ----------
    G : nx.Graph
        The bond graph in networkx format
    filename : str, optional
        The name of the file to save the bond graph. The default is 'nx_BG.html'.

    Returns
    -------
    None

    side effect
    ------------
    Save the bond graph in pyvis network format to html for visualization
    Use the default web browser to open the html file
        
    """

    pyvis_net = _pyvis_BG(G)
    _pyvis_net_fomat(pyvis_net)
    _nxBG_pyvis_net(pyvis_net, filename)
    browser = webbrowser.get()
    browser.open(os.path.abspath(filename))

def symExp_to_funcEval(symExp, df):
    """
    Evaluate the SymPy expression with the values in the dataframe

    Parameters
    ----------
    symExp : str
        The SymPy expression
    df : pandas.DataFrame
        The dataframe that contains the value of the SymPy expression.
        The column name of the dataframe should be the variable name in the SymPy expression.
        The column name of the dataframe should be 't' for the time variable.

    Returns
    -------
    func_eval : numpy.ndarray
        The numerical value of the SymPy expression
    func_eval_integrate : float
        The numerical value of the integral of the SymPy expression
    func_eval_integrate_abs : float
        The numerical value of the integral of the absolute value of the SymPy expression
    func_eval_integrate_cumulative : numpy.ndarray
        The numerical value of the cumulative integral of the SymPy expression
    """
    list_vars=list(symExp.free_symbols)
    list_vars_str=[str(var) for var in list_vars]
    func=lambdify(list_vars,symExp,'numpy')
    func_eval=func(*[df[var] for var in list_vars_str])
    func_eval_integrate=np.trapz(func_eval,df['t'])
    func_eval_integrate_abs=np.trapz(np.abs(func_eval),df['t'])
    func_eval_integrate_cumulative=integrate.cumulative_trapezoid(func_eval,x=df['t'],initial=0)

    return func_eval,func_eval_integrate,func_eval_integrate_abs,func_eval_integrate_cumulative

def calc_energy(P_comp_expr_dict,result_csv):
    """
    Calculate the energy and activity of the bond graph model

    Parameters
    ----------
    P_comp_expr_dict : str
        The dictionary of the expressions of the power of the components in the bond graph model
        The key of the dictionary is the component name
    result_csv : str
        The file path of the csv file to save the result

    Returns
    -------
    None

    side effect
    ------------
    Save the energy and activity of the bond graph model to a csv file

    """
    df_result_csv=pd.read_csv(result_csv)
    # get the path of the csv file
    csv_path=Path(result_csv).parent
    csv_file_name=Path(result_csv).stem
    csv_file_power_comp=csv_path/(csv_file_name+'_comp_power.csv')
    csv_file_activity_comp=csv_path/(csv_file_name+'_comp_activity.csv')
    E_comp_val=np.zeros(len(P_comp_expr_dict))
    A_comp_val=np.zeros(len(P_comp_expr_dict))
    AI_comp_val=np.zeros(len(P_comp_expr_dict))
    df_power_comp=pd.DataFrame(columns=P_comp_expr_dict.keys())
    df_activity_comp=pd.DataFrame(columns=P_comp_expr_dict.keys())
    df_power_comp['t']=df_result_csv['t']
    df_activity_comp['t']=df_result_csv['t']
    i=0
    A_total=0
    for comp in P_comp_expr_dict:
        P_comp_expr=P_comp_expr_dict[comp]
        df_power_comp[comp],E_comp_val[i],A_comp_val[i],df_activity_comp[comp]=symExp_to_funcEval(P_comp_expr,df_result_csv)
        A_total+=A_comp_val[i]
        i=i+1
        
    AI_comp_val=A_comp_val/A_total*100
    df_power_comp.to_csv(csv_file_power_comp,index=False)
    df_activity_comp.to_csv(csv_file_activity_comp,index=False)
    dict_activity={'Componentl list':list(P_comp_expr_dict.keys()),'Energy': E_comp_val.tolist(), 'Activity': A_comp_val.tolist(), 'Activity Index': AI_comp_val.tolist(),
                   }
    save_json(dict_activity,csv_path/(csv_file_name+'_activity.json'))



    
        
    