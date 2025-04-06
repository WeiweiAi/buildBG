from build_nxBG import getPowerPorts, getSignalPorts,load_nxBG_json,save_nxBG_json, nxBG_Energy
from sympy import *
import copy
import json
"""
This module is used to refine the bond graph model in the networkx format
The BondElement are refined based on the templates in the JSON file (./components/BG_components.json)
"""
def nxBG_refine_domain(G):
    """
    Specify the domain of the BondElement based on user inputs
    
    """
    BondElement_domain = ['electrical', 'mechanical', 'hydraulic', 'biochemical']
    # print the subclass with the numberings
    for idx, domain in enumerate(BondElement_domain, start=1):
        print(f"{idx}. {domain}")
    for node in G.nodes:
        if G.nodes[node]['a']=='BondElement' and 'domain' not in G.nodes[node].keys():
            print(f"The element: {node} belongs to domain:")
            subclass = input()
            G.nodes[node]['domain'] = BondElement_domain[int(subclass)-1]

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
        direction='out'
    elif any(data.get('relationship') == 'hasPowerPort' for u, v, data in G.in_edges(port,data=True)):
        direction='in'
    else:
        raise ValueError('The direction is not correct')
    return direction

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
        elif len(powerPorts)>=2:       
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
            if bg_comp_dict["ports"][powerPortN]['orientation']==_checkPortDirection(G,powerPort):
                pass
            elif bg_comp_dict["ports"][powerPortN]['orientation']=='in' and _checkPortDirection(G,powerPort)=='out':
                for u, v, data in G.out_edges(bondElement,data=True):
                    G.remove_edge(u, v)
                    G.add_edge(v,u,**data)
            elif bg_comp_dict["ports"][powerPortN]['orientation']=='out' and _checkPortDirection(G,powerPort)=='in':
                for u, v, data in G.in_edges(bondElement,data=True):
                    G.remove_edge(u, v)
                    G.add_edge(v,u,**data) 
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
        if G.nodes[node]['a']=='BondElement' and 'hasConstitutiveRelation' in G.nodes[node].keys():
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

if __name__ == "__main__": 
    
    json_file = './components/BG_components.json'
    bgComponents = json.load(open(json_file))
    nx_BG_file = './data/nx_BG_domain.json'
    G = load_nxBG_json(nx_BG_file)
    nxBG_refine_components(G,bgComponents)
    nx_BG_file = './data/nx_BG_refine.json'
    save_nxBG_json(G, nx_BG_file)
    nxBG_Energy(G)
    nx_BG_file = './data/nx_BG_energy.json'
    save_nxBG_json(G, nx_BG_file)
    voi={'propertyName': 't', 'units': 'second'}
    nxBG_refine_constitutive_relations(G)   
    nx_BG_file = './data/nx_BG_constitutive_energy.json'
    save_nxBG_json(G, nx_BG_file)
