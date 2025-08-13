from build_nxBG import getPowerPorts,load_nxBG_json,save_nxBG_json, nxBG_Energy,save_nxBG_html
from sympy import *
import copy
import json
import pandas as pd
"""
This module is used to refine the bond graph model in the networkx format
The BondElement are refined based on the templates in the JSON file (./components/BG_components.json)
"""

def nxBG_refine_domain(G,BG_domain):
    """
    Refine the domain of the BondElement and multiJunction in the bond graph 
    
    """
    for node in G.nodes:
        if G.nodes[node]['a']=='BondElement' or (G.nodes[node]['a']=='JunctionStructure' and (G.nodes[node]['subClass']=='TF' or G.nodes[node]['subClass']=='GY')):
           G.nodes[node]['domain']=BG_domain[G.nodes[node]['domain']]

def _checkPortDirection(G,port):
    """
    Check the direction of the port, either 'in' or 'out'
    'in' means the power is positive when it enters the bond element from the port
    'out' means the power is positive when it leaves the bond element to the port

    parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx
    port : str
        The name of the port

    returns
    -------
    direction : str
        The direction of the port, either 'in' or 'out'

    """
    if any(data.get('relationship') == 'hasPowerPort' for u, v, data in G.out_edges(port,data=True)):
        direction='in'
    elif any(data.get('relationship') == 'hasPowerPort' for u, v, data in G.in_edges(port,data=True)):
        direction='out'
    else:
        raise ValueError('The direction is not correct')
    return direction

def _flipEdge(G, port, direction):
    """
    Flip the direction of the edge connected to the port
    """
    if direction == 'in':  # from in to out
        edges_to_flip = [(u, v, data) for u, v, data in G.in_edges(port, data=True) if data.get('a') == 'PowerBond']
        for u, v, data in edges_to_flip:
            G.add_edge(v, u, **data)
            G.remove_edge(u, v)
    elif direction == 'out':  # from out to in
        edges_to_flip = [(u, v, data) for u, v, data in G.out_edges(port, data=True) if data.get('a') == 'PowerBond']
        for u, v, data in edges_to_flip:
            G.add_edge(v, u, **data)
            G.remove_edge(u, v)


def nxBG_refine_component(G,bondElement,bgComponents):
    """
    Refine the components in the bond graph

    parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx
    bondElement : str
        The name of the bond element
    bgComponents : dict
        The dictionary of the bond graph components (template)

    returns
    -------
    none

    side effects
    ------------
    The modelParameter and modelState of the bond element are updated
    The power ports are updated with the propertyName and symbol
    The constitutive relations are added if the bond element has any
    The orientation of the power ports are checked and updated if necessary
    """
    def check_sharedVar(G,bondElement,powerPorts, causality): 
        if len(powerPorts)==1:
            return True
        elif len(powerPorts)>=2:    # e.g., f_0==f_1 in Re component
            output_var_0=f'{causality[0]}_{powerPorts[0].split("_")[-1]}'
            output_var_1=f'{causality[0]}_{powerPorts[1].split("_")[-1]}'
            for expr in G.nodes[bondElement]['hasConstitutiveRelation']:
                var_=solve(expr, symbols(output_var_0))
                if len(var_)==1 and var_[0]==symbols(output_var_1):
                    return True
            return False   
           
    powerPorts=getPowerPorts(G,bondElement)
    bg_domain_dict={}
    bg_comp_dict={}
    for key in bgComponents:
        if bgComponents[key]['domain'] == G.nodes[bondElement]['domain']:
            bg_domain_dict=bgComponents[key]
            break    
    if len(bg_domain_dict)==0:
        raise ValueError('The domain is unspecified!')
    else:
        for component in bg_domain_dict['components'].values():
            if G.nodes[bondElement]['subClass']==component['metamodel']:
                bg_comp_dict=copy.deepcopy(component)
                break
    if len(bg_comp_dict)==0:
        raise ValueError('The component is not found!')
    if 'params' in bg_comp_dict:
        G.nodes[bondElement]['modelParameter']=copy.deepcopy(bg_comp_dict['params'])
        for param in bg_comp_dict['params']:
            quantity=G.nodes[bondElement]['modelParameter'][param]
            if ('physical constants' in bg_domain_dict and param in bg_domain_dict['physical constants'].keys()) or ('thermodynamic parameters' in bg_domain_dict and param in bg_domain_dict['thermodynamic parameters'].keys()):               
                quantity['propertyName']=quantity['symbol']
            else:
                quantity['propertyName']=quantity['symbol']+'_'+bondElement
    if 'constitutive_relations' in bg_comp_dict:
        G.nodes[bondElement]['hasConstitutiveRelation']=[]
        G.nodes[bondElement]['hasConstitutiveRelation']=copy.deepcopy(bg_comp_dict['constitutive_relations'])
    if  G.nodes[bondElement]['subClass']=='C' or G.nodes[bondElement]['subClass']=='MC' or G.nodes[bondElement]['subClass']=='E' or G.nodes[bondElement]['subClass']=='ME':
        G.nodes[bondElement]['modelState']={}
        G.nodes[bondElement]['modelState']['q_0']=copy.deepcopy(bg_domain_dict['generalized displacement'])
        quantity=G.nodes[bondElement]['modelState']['q_0']
        quantity['propertyName']=quantity['symbol']+'_'+bondElement
    if  G.nodes[bondElement]['subClass']=='I' or G.nodes[bondElement]['subClass']=='MI':
        G.nodes[bondElement]['modelState']={}
        G.nodes[bondElement]['modelState']['p_0']=copy.deepcopy(bg_domain_dict['generalized momentum'])
        quantity=G.nodes[bondElement]['modelState']['p_0']
        quantity['propertyName']=quantity['symbol']+'_'+bondElement
    if  G.nodes[bondElement]['subClass']=='IC' or G.nodes[bondElement]['subClass']=='MIC' or G.nodes[bondElement]['subClass']=='IE' or G.nodes[bondElement]['subClass']=='MIE':
        G.nodes[bondElement]['modelState']={}
        G.nodes[bondElement]['modelState']['q_0']=copy.deepcopy(bg_domain_dict['generalized displacement'])
        quantity=G.nodes[bondElement]['modelState']['q_0']
        quantity['propertyName']=quantity['symbol']+'_'+bondElement
        G.nodes[bondElement]['modelState']['p_1']=copy.deepcopy(bg_domain_dict['generalized momentum'])           
        quantity2=G.nodes[bondElement]['modelState']['p_1']
        quantity2['propertyName']=quantity2['symbol']+'_'+bondElement
    if len(powerPorts)!=len(bg_comp_dict["ports"]):
        raise ValueError('The number of ports does not match.') 
    else:
        for powerPort in powerPorts:
            powerPortN=powerPort.split('_')[-1] # The format of the power port is 'bondElement_powerPortN'
            G.nodes[powerPort]['effort']=copy.deepcopy(bg_domain_dict['effort'])
            G.nodes[powerPort]['flow']=copy.deepcopy(bg_domain_dict['flow'])
            effort=G.nodes[powerPort]['effort']
            flow=G.nodes[powerPort]['flow']
            if check_sharedVar(G,bondElement,powerPorts,'effort'):
                effort['propertyName']=effort['symbol']+'_'+bondElement
            else:
                effort['propertyName']=effort['symbol']+'_'+powerPort
            if check_sharedVar(G,bondElement,powerPorts,'flow'):
                flow['propertyName']=flow['symbol']+'_'+bondElement
            else:
                flow['propertyName']=flow['symbol']+'_'+powerPort
            # check the orientation of the power port           
            if bg_comp_dict["ports"][powerPortN]['orientation']=='in':
                _flipEdge(G,powerPort,'out')# align the orientation of the bond and the power port
                if _checkPortDirection(G,powerPort)=='out':
                    data=G.get_edge_data(bondElement,powerPort)
                    G.remove_edge(bondElement,powerPort)
                    G.add_edge(powerPort,bondElement,**data)
            elif bg_comp_dict["ports"][powerPortN]['orientation']=='out':
                _flipEdge(G,powerPort,'in')# align the orientation of the bond and the power port
                if _checkPortDirection(G,powerPort)=='in':
                    data=G.get_edge_data(powerPort,bondElement)
                    G.remove_edge(powerPort,bondElement)
                    G.add_edge(bondElement,powerPort,**data)
            else:
                 raise ValueError('The orientation does not match.')    

def nxBG_refine_multiJunc(G,multiJunc,bgComponents):
    """
    Refine the components in the bond graph

    parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx
    multiJunc : str
        The name of the junction
    bgComponents : dict
        The dictionary of the bond graph components (template)

    returns
    -------
    none

    side effects
    ------------
    The modelParameter the junction are updated
    The power ports are updated with the propertyName and symbol
    The constitutive relations are added if the junction has any

    """
    powerPorts=getPowerPorts(G,multiJunc)
    bg_domain_dict={}
    bg_comp_dict={}
    for key in bgComponents:
        if bgComponents[key]['domain'] == G.nodes[multiJunc]['domain']:
            bg_domain_dict=bgComponents[key]
            break    
    if len(bg_domain_dict)==0:
        raise ValueError('The domain is unspecified!')
    else:
        for component in bg_domain_dict['components'].values():
            if G.nodes[multiJunc]['subClass']==component['metamodel']:
                bg_comp_dict=copy.deepcopy(component)
                break
    if len(bg_comp_dict)==0:
        raise ValueError('The component is not found!')
    if 'params' in bg_comp_dict and 'modelParameter' not in G.nodes[multiJunc].keys():
        G.nodes[multiJunc]['modelParameter']={}
        G.nodes[multiJunc]['modelParameter']=copy.deepcopy(bg_comp_dict['params'])
        for param in bg_comp_dict['params']:
            quantity=G.nodes[multiJunc]['modelParameter'][param]
            if ('physical constants' in bg_domain_dict and param in bg_domain_dict['physical constants'].keys()) or ('thermodynamic parameters' in bg_domain_dict and param in bg_domain_dict['thermodynamic parameters'].keys()):               
                quantity['propertyName']=quantity['symbol']
            else:
                quantity['propertyName']=quantity['symbol']+'_'+multiJunc
    if len(powerPorts)!=len(bg_comp_dict["ports"]):
        raise ValueError('The number of ports does not match.') 
    else:
        for powerPort in powerPorts:
            powerPortN=powerPort.split('_')[-1]          
             # check the orientation of the power port
            if bg_comp_dict["ports"][powerPortN]['orientation']=='in':
                if _checkPortDirection(G,powerPort)=='out':
                    data=G.get_edge_data(multiJunc,powerPort)
                    G.remove_edge(multiJunc,powerPort)
                    G.add_edge(powerPort,multiJunc,**data)
            elif bg_comp_dict["ports"][powerPortN]['orientation']=='out':
                if _checkPortDirection(G,powerPort)=='in':
                    data=G.get_edge_data(powerPort,multiJunc)
                    G.remove_edge(powerPort,multiJunc)
                    G.add_edge(multiJunc,powerPort,**data)
            else:
                 raise ValueError('The orientation does not match.')

            
def nxBG_refine_components(G,bgComponents):
    
    """
    This function is used to refine the components in the bond graph

    Parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx
        Nodes: BondElement, JunctionStructure, PowerPort, SignalPort
        Edges: PowerBond, SignalBond, hasPowerPort, hasSignalPort
    
    bgComponents : dict
        The dictionary of the bond graph components (template)
    
    """
    for node in G.nodes:
        if G.nodes[node]['a']=='BondElement':
            nxBG_refine_component(G,node,bgComponents)
        elif G.nodes[node]['a']=='JunctionStructure' and (G.nodes[node]['subClass']=='TF' or G.nodes[node]['subClass']=='GY'):
            nxBG_refine_multiJunc(G,node,bgComponents)

def update_expr(expr,var_dict):
    for ikey,ivar in var_dict.items():
        if 'expression' in var_dict[ikey].keys():
            expr=expr.subs(symbols(ikey),sympify(var_dict[ikey]['expression']))
        else:
            expr=expr.subs(symbols(ikey),symbols(var_dict[ikey]['propertyName']))
    return expr

def nxBG_refine_constitutive_relations(G):

    """
    Refine the constitutive relations in the bond graph
    The constitutive relations are in the form of a list of strings
    The effort and flow variables are in the form of 'e_0', 'f_0', 'e_1', 'f_1', etc.
    '0' and '1' are the port numbers

    Parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx       
    """

    for node in G.nodes: 
        if 'hasConstitutiveRelation' in G.nodes[node].keys():
            G.nodes[node]['constitutive_eqs']={}
            var_dict={}
            powerPorts=getPowerPorts(G,node)
            if 'modelParameter' in G.nodes[node].keys():
                var_dict.update(copy.deepcopy(G.nodes[node]['modelParameter']))
            if 'modelState' in G.nodes[node].keys():
                var_dict.update (copy.deepcopy(G.nodes[node]['modelState']))
            for powerPort in powerPorts:
                powerPortN=powerPort.split('_')[-1]
                var_dict.update({f'e_{powerPortN}':G.nodes[powerPort]['effort']})
                var_dict.update({f'f_{powerPortN}':G.nodes[powerPort]['flow']})
            for str_expr in G.nodes[node]['hasConstitutiveRelation']:
                if 'ode' in str_expr:
                    yvar=str_expr.split('ode(')[1].split(',')[0] # the form of the ode is 'ode(yvar, t)-expr'
                    other_terms=str_expr.split('ode(')[1].split(')')[1]
                    expr=-sympify(other_terms) # assume the ode is on the left hand side
                    expr=update_expr(expr,var_dict)
                    LHS_name=var_dict[yvar]['propertyName']
                    if LHS_name not in G.nodes[node]['constitutive_eqs'].keys():
                        G.nodes[node]['constitutive_eqs'][LHS_name]=(ccode(expr), 't')
                else:
                    expr=sympify(str_expr)
                    for powerPort in powerPorts:
                        powerPortN=powerPort.split('_')[-1]
                        causality=G.nodes[powerPort]['causality']
                        output_var=f'{causality[0]}_{powerPortN}'
                        if symbols(output_var) in expr.free_symbols:
                            var_=solve(expr, symbols(output_var))
                            if len(var_)==0:
                                raise ValueError('The equation is not solvable')
                            else:
                                solved_expr=var_[0]
                                solved_expr=update_expr(solved_expr,var_dict)
                                LHS_name=var_dict[output_var]['propertyName']
                                if LHS_name not in G.nodes[node]['constitutive_eqs'].keys():
                                    G.nodes[node]['constitutive_eqs'][LHS_name]=(ccode(solved_expr), '')

def nxBG_getParameters(G):
    """
    Get the parameters from the bond graph

    Parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx

    Returns
    -------
    parameters : dict
        The parameters in the bond graph

    """
    parameters={}
    for node in G.nodes:
        if 'modelParameter' in G.nodes[node].keys():
            for param in G.nodes[node]['modelParameter']:
                if 'value' in G.nodes[node]['modelParameter'][param].keys():
                    parameters[G.nodes[node]['modelParameter'][param]['propertyName']]=G.nodes[node]['modelParameter'][param]['value']
                else:
                    parameters[G.nodes[node]['modelParameter'][param]['propertyName']]=1.0 # default value
        if 'modelState' in G.nodes[node].keys():
            for param in G.nodes[node]['modelState']:
                if 'value' in G.nodes[node]['modelState'][param].keys():
                    parameters[G.nodes[node]['modelState'][param]['propertyName']]=G.nodes[node]['modelState'][param]['value']
                else:
                    parameters[G.nodes[node]['modelState'][param]['propertyName']]=1.0 # default value
    return parameters

def nxBG_refine_parameters(G,parameters):
    """
    Refine the parameters in the bond graph

    Parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx
    parameters : dict
        The parameters in the bond graph

    Returns
    -------
    none

    side effects
    ------------
    The parameters in the bond graph are updated with the values from the dictionary

    """
    
    for node in G.nodes:
        if 'modelParameter' in G.nodes[node].keys():
            for param in G.nodes[node]['modelParameter']:
                if G.nodes[node]['modelParameter'][param]['propertyName'] in parameters.keys():
                    G.nodes[node]['modelParameter'][param]['value']=parameters[G.nodes[node]['modelParameter'][param]['propertyName']]
        if 'modelState' in G.nodes[node].keys():
            for param in G.nodes[node]['modelState']:
                if G.nodes[node]['modelState'][param]['propertyName'] in parameters.keys():
                    G.nodes[node]['modelState'][param]['value']=parameters[G.nodes[node]['modelState'][param]['propertyName']]

def nxBG_refine_paramCsv(G,bg_params_csv):
    """
    Refine the parameters in the bond graph

    Parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx
    bg_params_csv : str
        The path of the csv file containing the parameters

    Returns
    -------
    none

    side effects
    ------------
    The parameters in the bond graph are updated with the values from the csv file

    """
    
    parameters=pd.read_csv(bg_params_csv, header=None).set_index(0).T.to_dict('records')[0]
    nxBG_refine_parameters(G,parameters)

def nxBG_getParameters_csv(G,bg_params_csv):
    """
    Get the parameters from the bond graph and write them to a csv file

    Parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx
    bg_params_csv : str
        The path of the csv file to save the parameters

    Returns
    -------
    none

    side effects
    ------------
    The parameters in the bond graph are written to the csv file

    """
    
    parameters=nxBG_getParameters(G)
    pd.DataFrame.from_dict(parameters, orient='index').to_csv(bg_params_csv, header=False)

def nxBG_refine(nxBGJson,bgComponentsJson,BG_domainJson,nxBGRefinedJson=None,bg_params_csv=None):
    """
    Refine the bond graph model in the networkx format
    The BondElement are refined based on the templates in the JSON file (./components/BG_components.json)

    """
    
    G = load_nxBG_json(nxBGJson)
    bgComponents = json.load(open(bgComponentsJson))
    BG_domain = json.load(open(BG_domainJson))
    
    nxBG_refine_domain(G,BG_domain)    
    nxBG_refine_components(G,bgComponents)
    nxBG_Energy(G)
    nxBG_refine_constitutive_relations(G)
    if bg_params_csv:
        nxBG_refine_paramCsv(G,bg_params_csv)
    if nxBGRefinedJson:
        save_nxBG_json(G, nxBGRefinedJson)
        save_nxBG_html(G, nxBGRefinedJson.split('.json')[0]+'.html')
    else:
        nx_BG_refined_file = nxBGJson.split('.json')[0]+'_refined.json'
        save_nxBG_json(G, nx_BG_refined_file)
        save_nxBG_html(G, nx_BG_refined_file.split('.json')[0]+'.html')

if __name__ == "__main__": 
    
    bgComponents = './components/BG_components.json'
    BG_domain_file = './components/BG_domain.json'
    nx_BG_file = './data/nx_BG.json'
    
    nxBG_refine(nx_BG_file,bgComponents,BG_domain_file,bg_params_csv='./data/bg_params.csv')
    nx_BG_file = nx_BG_file.split('.json')[0]+'_refined.json'
    G=load_nxBG_json(nx_BG_file)
    path_='./data/'
    save_nxBG_html(G, path_+'nx_BG_refined.html')



