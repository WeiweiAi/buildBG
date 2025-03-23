import networkx as nx
from buildBG import load_matrix
import matplotlib.pyplot as plt
from networkx.readwrite import json_graph
import json
from pyvis.network import Network
from sympy import *

def nxBG(G, matrix, direction='e2f', connection='PowerBond'):
    """
    This function is used to build the bond graph using networkx based on stoichiometric matrix
    Note: The 0 or 1 junctions may not be added to the bond graph, which can be added by nxBG_addJunction function later
    
    Parameters
    ----------
    G : nx.DiGraph
        The graph to build the bond graph
    matrix : str
        The path of the matrix file
    direction : str
        The direction of the bond graph, either 'e2f' or 'f2e'
    connection : str
        The connection type, either 'PowerBond' or 'SignalBond'

    Returns
    -------
    G : nx.DiGraph
        The bond graph built using networkx
        Nodes: BondElement, JunctionStructure, PowerPort, SignalPort
        Edges: PowerBond, SignalBond, hasPowerPort, hasSignalPort
    """
    eName, eID, ePort, fName, fID, fPort,N=load_matrix(matrix)
    if direction=='f2e':
        N=-N
    cName=eName+fName # component name, should be unique for the graph
    cPort=ePort+fPort # port name, unique for each component
    cID=eID+fID # can be used to specify the component type, such as resistor, capacitor, etc. '1' and '0' are used for the junctions
    for i in range(len(cName)):
        if cName[i] not in G.nodes:
            if cID[i]=='1': # BondGraph Junction
                G.add_node(cName[i], a='JunctionStructure', subClass='OneJunctionStructure')
            elif cID[i]=='0': # BondGraph Junction
                G.add_node(cName[i], a='JunctionStructure', subClass='ZeroJunctionStructure')
            elif cID[i]=='TF': # BondGraph Transformer
                G.add_node(cName[i], a='JunctionStructure', subClass='Transformer')
            elif cID[i]=='GY': # BondGraph Gyrator
                G.add_node(cName[i], a='JunctionStructure', subClass='Gyrator')
            else: # BondGraph Element
                G.add_node(cName[i], a='BondElement',component=cID[i])
        iport=f'{cName[i]}_{cPort[i]}'
        if connection=='PowerBond':
            if i<len(eName):
                G.add_node(iport, a='PowerPort', isPortOf=G.nodes[cName[i]]['a'], effort=None, flow=None, causality='effort')
            if i>=len(eName):
                G.add_node(iport, a='PowerPort', isPortOf=G.nodes[cName[i]]['a'], effort=None, flow=None, causality='flow')
            
            G.add_edge(iport, cName[i], relationship='hasPowerPort')# By default, the positive reference direction for energy flows is inward
        elif connection=='SignalBond':
            G.add_node(iport, a='SignalPort', isPortOf=G.nodes[cName[i]]['a'], signal=None)
            G.add_edge(iport, cName[i], relationship='hasSignalPort')
        else:
            raise ValueError('The connection type is not correct')
        
    for i in range(len(eName)):
        for j in range(len(fName)):
            iport=f'{eName[i]}_{ePort[i]}'
            jport=f'{fName[j]}_{fPort[j]}'
            if N[i][j]>0:
                G.add_edge(iport,jport, multiplier=str(N[i][j]), a=connection, effort=None, flow=None)
            elif N[i][j]<0:   
                G.add_edge(jport,iport, multiplier=str(-N[i][j]), a=connection, effort=None, flow=None)
    return G

def nxBG_addJunction(G):
    """
    This function is used to add the 0 or 1 junctions to the bond graph 
    in case the PowerPort of the BondElement have more than 2 powerbond connections
    
    Parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx
        Nodes: BondElement, JunctionStructure, PowerPort, SignalPort
        Edges: PowerBond, SignalBond, hasPowerPort, hasSignalPort
        
    Returns
    -------
    None
    """
    PowerPort_BondElement=[node for node in G.nodes if G.nodes[node]['a']=='PowerPort' and G.nodes[node]['isPortOf']=='BondElement']
    for node in PowerPort_BondElement:        
        if G.degree(node)>2:
            jName='J_'+node
            G.add_node(jName, a='JunctionStructure')
            jPort=0
            # move the edges to the junction
            out_edges=[(u, v) for u, v, data in G.out_edges(node,data=True) if data.get('a') == 'PowerBond']
            in_edges=[(u, v) for u, v, data in G.in_edges(node,data=True) if data.get('a') == 'PowerBond']
            for edge in out_edges:
                jPort+=1
                jport=f'{jName}_{jPort}'
                G.add_node(jport, a='PowerPort', isPortOf=G.nodes[jName]['a'], effort=None, flow=None, causality=G.nodes[node]['causality'])
                G.add_edge(jName,jport, relationship='hasPowerPort')
                G.add_edge(jport,edge[1], multiplier=G.edges[edge]['multiplier'], a=G.edges[edge]['a'], effort=None, flow=None)
                G.remove_edge(edge[0],edge[1])
            for edge in in_edges:
                jPort+=1
                jport=f'{jName}_{jPort}'
                G.add_node(jport, a='PowerPort', isPortOf=G.nodes[jName]['a'], effort=None, flow=None, causality=G.nodes[node]['causality'])
                G.add_edge(jport,jName, relationship='hasPowerPort')
                G.add_edge(edge[0],jport, multiplier=G.edges[edge]['multiplier'], a=G.edges[edge]['a'], effort=None, flow=None)
                G.remove_edge(edge[0],edge[1]) 
            jPort=0
            jport=f'{jName}_{jPort}'
            if G.nodes[node]['causality']=='effort':
                G.nodes[jName]['subClass']='ZeroJunctionStructure'
                G.add_node(jport, a='PowerPort', isPortOf=G.nodes[jName]['a'], effort=None, flow=None, causality='flow')
                G.add_edge(jport,node, multiplier='1', a='PowerBond', effort=None, flow=None)                
            elif G.nodes[node]['causality']=='flow':
                G.nodes[jName]['subClass']='OneJunctionStructure'
                G.add_node(jport, a='PowerPort', isPortOf=G.nodes[jName]['a'], effort=None, flow=None, causality='effort')
                G.add_edge(jport,node,multiplier='1', a='PowerBond', effort=None, flow=None)
            else:
                raise ValueError('The causality is not correct.')      
            G.add_edge(jName,jport, relationship='hasPowerPort')

def nxBG_addJunction_multi(G):
    """
    This function is used to replace the PowerBond with the Transformer in case the PowerBond has the multiplier not equal to 1

    Parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx
        Nodes: BondElement, JunctionStructure, PowerPort, SignalPort
        Edges: PowerBond, SignalBond, hasPowerPort, hasSignalPort

    Returns
    -------
    None
    
    """
    PowerBonds=[edge for edge in G.edges if 'a' in G.edges[edge].keys() and  G.edges[edge]['a']=='PowerBond' and G.edges[edge]['multiplier']!='1']    
    for edge in PowerBonds:
        source=edge[0]
        target=edge[1]
        TFName='TF_'+source+'_'+target
        G.add_node(TFName, a='JunctionStructure', subClass='Transformer', 
                   modelParameter={'param_1':
                                    {'propertyValue':G.edges[edge]['multiplier'], 'units':'dimensionless',
                                    'propertyName':TFName}
                                    },           
            )
        # add the power ports of the transformer
        port_0=f'{TFName}_0'
        port_1=f'{TFName}_1'
        G.add_node(port_0, a='PowerPort', isPortOf=G.nodes[TFName]['a'], effort=None, flow=None, causality='flow')
        G.add_edge(TFName,port_0, relationship='hasPowerPort')
        G.add_node(port_1, a='PowerPort', isPortOf=G.nodes[TFName]['a'], effort=None, flow=None, causality='effort')
        G.add_edge(TFName,port_1, relationship='hasPowerPort')
        if G.nodes[source]['causality']=='effort':
            G.add_edge(source,port_0, multiplier='1', a='PowerBond', effort=None, flow=None)
            G.add_edge(port_1,target, multiplier='1', a='PowerBond', effort=None, flow=None)
        elif G.nodes[source]['causality']=='flow':
            G.add_edge(source,port_1, multiplier='1', a='PowerBond',  effort=None, flow=None)
            G.add_edge(port_0,target, multiplier='1', a='PowerBond',  effort=None, flow=None)
        else:
            raise ValueError('The causality is not correct')
        G.remove_edge(source,target)

def nxBG_addJunctions(G):
    """
    This function is used to add transformers, 0 and 1 junctions
    """
    nxBG_addJunction_multi(G)
    nxBG_addJunction(G)

def getPowerPorts(G,bondElement):
    """
    This function is used to get the power ports of the bond element

    Parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx
        Nodes: BondElement, JunctionStructure, PowerPort, SignalPort
        Edges: PowerBond, SignalBond, hasPowerPort, hasSignalPort
    bondElement : str
        The name of the bond element

    Returns
    -------
    ports : list
        The list of the power ports of the bond element
    """
    ports=[out_edge[1] for out_edge in G.out_edges(bondElement) if 'relationship' in G.edges[out_edge] and G.edges[out_edge]['relationship']=='hasPowerPort'] # get the power ports of the bond element
    ports+= [in_edge[0] for in_edge in G.in_edges(bondElement) if 'relationship' in G.edges[in_edge] and G.edges[in_edge]['relationship']=='hasPowerPort']
    return ports

def getSignalPorts(G,bondElement):
    """
    This function is used to get the signal ports of the bond element

    Parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx
        Nodes: BondElement, JunctionStructure, PowerPort, SignalPort
        Edges: PowerBond, SignalBond, hasPowerPort, hasSignalPort
    bondElement : str
        The name of the bond element

    Returns
    -------
    ports : list
        The list of the signal ports of the bond element
    """
    ports=[out_edge[1] for out_edge in G.out_edges(bondElement) if 'relationship' in G.edges[out_edge] and G.edges[out_edge]['relationship']=='hasSignalPort'] # get the signal ports of the bond element
    return ports

def nxBG_initEnergy(G):
    """
    This function is used to initialize:
    the output variable of the PowerPort of the BondElement based on the causality of the PowerBond, either effort or flow,
    the variable is stored as 'propertyName'

    """
    for node in G.nodes:
        if G.nodes[node]['a']=='PowerPort' and G.nodes[node]['isPortOf']=='BondElement':
            if G.nodes[node]['causality']=='effort' and G.nodes[node]['effort'] is None:
                G.nodes[node]['effort']={}
                G.nodes[node]['effort']['propertyName']='e_'+ node
            if G.nodes[node]['causality']=='flow' and G.nodes[node]['flow'] is None:
                G.nodes[node]['flow']={}
                G.nodes[node]['flow']['propertyName']='f_'+ node

def viaPowerbond(G,edge):
    """
    This function is used to propagate the energy through the PowerBond
    
    Parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx
        Nodes: BondElement, JunctionStructure, PowerPort, SignalPort
        Edges: PowerBond, SignalBond, hasPowerPort, hasSignalPort
    edge : tuple
        The edge of the PowerBond

    Returns
    -------
    Boolean
        True if the energy is propagated through the PowerBond, otherwise False

    side effects
    ------------
    The expression of the effort and flow of the PowerPorts and the PowerBond are updated if the energy is propagated through the PowerBond

    Raises
    ------
    ValueError
        If the causality of the PowerPort is not correct
        If the expression of the PowerPort is not consistent with the PowerBond
    
    """    
    source=edge[0]
    target=edge[1]
    source_var_out=G.nodes[source]['causality']
    target_var_out=G.nodes[target]['causality']
    if source_var_out==target_var_out:
        raise ValueError('The causality is not correct')
    else:
        if G.nodes[source][source_var_out] is not None:
            if G.edges[edge][source_var_out] is None:
                G.edges[edge][source_var_out]={}
            if 'propertyName' in G.nodes[source][source_var_out].keys(): 
                if 'expression' in G.edges[edge][source_var_out].keys():
                    if symbols(G.nodes[source][source_var_out]['propertyName'])!=G.edges[edge][source_var_out]['expression']:
                        raise ValueError(f'The {source_var_out} of the port {source} is not consistent with the PowerBond')
                else:                   
                    G.edges[edge][source_var_out]['expression']=symbols(G.nodes[source][source_var_out]['propertyName'])
            elif 'expression' in G.nodes[source][source_var_out].keys():
                if 'expression' in G.edges[edge][source_var_out].keys():
                    if G.nodes[source][source_var_out]['expression']!=G.edges[edge][source_var_out]['expression']:
                        raise ValueError(f'The {source_var_out} of the port {source} is not consistent with the PowerBond')
                else:
                    G.edges[edge][source_var_out]['expression']=G.nodes[source][source_var_out]['expression']
            else:
                raise ValueError(f'The {source_var_out} of the port {source} is unknown.')
            if G.nodes[target][source_var_out] is None: 
                G.nodes[target][source_var_out]={}
            if 'expression' in G.nodes[target][source_var_out].keys():
                if G.nodes[target][source_var_out]['expression']!=G.edges[edge][source_var_out]['expression']:
                    raise ValueError(f'The {source_var_out} of the port {target} is not consistent with the PowerBond')
            else:
                G.nodes[target][source_var_out]['expression']=G.edges[edge][source_var_out]['expression']       
        if G.nodes[target][target_var_out] is not None:
            if G.edges[edge][target_var_out] is None:
                G.edges[edge][target_var_out]={}
            if 'propertyName' in G.nodes[target][target_var_out].keys():
                if 'expression' in G.edges[edge][target_var_out].keys():
                    if symbols(G.nodes[target][target_var_out]['propertyName'])!=G.edges[edge][target_var_out]['expression']:
                        raise ValueError(f'The {target_var_out} of the port {target} is not consistent with the PowerBond')
                else:                    
                    G.edges[edge][target_var_out]['expression']=symbols(G.nodes[target][target_var_out]['propertyName'])
            elif 'expression' in G.nodes[target][target_var_out].keys():
                if 'expression' in G.edges[edge][target_var_out].keys():
                    if G.nodes[target][target_var_out]['expression']!=G.edges[edge][target_var_out]['expression']:
                        raise ValueError(f'The {target_var_out} of the port {target} is not consistent with the PowerBond')
                else:
                    G.edges[edge][target_var_out]['expression']=G.nodes[target][target_var_out]['expression']
            else:
                raise ValueError(f'The {target_var_out} of the port {target} is unknown.')
            if G.nodes[source][target_var_out] is None:
                G.nodes[source][target_var_out]={}
            if 'expression' in G.nodes[source][target_var_out].keys():
                if G.nodes[source][target_var_out]['expression']!=G.edges[edge][target_var_out]['expression']:
                    raise ValueError(f'The {target_var_out} of the port {source} is not consistent with the PowerBond')
            else:
                G.nodes[source][target_var_out]['expression']=G.edges[edge][target_var_out]['expression']

    if G.nodes[source][source_var_out] is None and G.nodes[target][target_var_out] is None:
        return False
    else:
        return True

def viaJunction(G,nodeJunction):
    """
    This function is used to propagate the energy through the 0 or 1 Junction

    Parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx
        Nodes: BondElement, JunctionStructure, PowerPort, SignalPort
        Edges: PowerBond, SignalBond, hasPowerPort, hasSignalPort
    nodeJunction : str
        The name of the junction
    """
    ports=getPowerPorts(G,nodeJunction) # get the power ports of the junction
    if G.nodes[nodeJunction]['subClass']=='ZeroJunctionStructure':
        shared_var_type='effort' # one effort in, multiple effort out
        compute_var_type='flow' # multiple flow in, one flow out
    elif G.nodes[nodeJunction]['subClass']=='OneJunctionStructure':
        shared_var_type='flow' # one flow in, multiple flow out
        compute_var_type='effort' # multiple effort in, one effort out  
    else:
        raise ValueError('The subclass is not correct') 
    
    for port in ports:
        if G.nodes[port][shared_var_type] is not None and G.nodes[port]['causality']==compute_var_type: # input to the junction
            if 'expression' in G.nodes[port][shared_var_type].keys():
                shared_var=G.nodes[port][shared_var_type]['expression']
                compute_direction=checkBondDirection(G,port)
                compute_port=port                                    
                break
            else:
               raise ValueError(f'The {shared_var_type} of the port {port} is unknown.')
        elif G.nodes[port][shared_var_type] is None and G.nodes[port]['causality']==compute_var_type: # input to the junction not solved
            return False
        
    compute_var_init=0
    known_compute_var=0
    for  port in ports:
        if G.nodes[port]['causality']==shared_var_type:
            if G.nodes[port][shared_var_type] is None:
                G.nodes[port][shared_var_type]={}
            if 'expression' in G.nodes[port][shared_var_type].keys():
                if G.nodes[port][shared_var_type]['expression']!=shared_var:
                    raise ValueError(f'The {shared_var_type} of the port {port} is not consistent with the junction')    
            else:
                G.nodes[port][shared_var_type]['expression']=shared_var
            # get the powerbond of the port and propagate the energy
            powerBonds= [u for u in G.out_edges(port) if 'a' in G.edges[u].keys() and G.edges[u]['a']=='PowerBond']+ [u for u in G.in_edges(port) if 'a' in G.edges[u].keys() and G.edges[u]['a']=='PowerBond']
            viaPowerbond(G,powerBonds[0])
            if G.nodes[port][compute_var_type] is not None:
                if 'expression' in G.nodes[port][compute_var_type].keys():
                    known_compute_var+=1
                    if checkBondDirection(G,port)==compute_direction: 
                        compute_var_init-=G.nodes[port][compute_var_type]['expression']
                    else:
                        compute_var_init+=G.nodes[port][compute_var_type]['expression']
                else:
                    raise ValueError(f'The {compute_var_type} of the port {port} is unknown.')
    
    if known_compute_var==len(ports)-1:
        if G.nodes[compute_port][compute_var_type] is None:
            G.nodes[compute_port][compute_var_type]={}
        if 'expression' in G.nodes[compute_port][compute_var_type].keys():
            if G.nodes[compute_port][compute_var_type]['expression']!=compute_var_init:
                raise ValueError(f'The {compute_var_type} of the port {compute_port} is not consistent with the junction')
        else:
            G.nodes[compute_port][compute_var_type]['expression']=compute_var_init
        powerBonds= [u for u in G.out_edges(compute_port) if 'a' in G.edges[u].keys() and G.edges[u]['a']=='PowerBond']+ [u for u in G.in_edges(compute_port) if 'a' in G.edges[u].keys() and G.edges[u]['a']=='PowerBond']
        viaPowerbond(G,powerBonds[0])
        return True
    else:
        return False

def viaTransformer(G,TF):
    """
    This function is used to propagate the energy through the Transformer

    Parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx
        Nodes: BondElement, JunctionStructure, PowerPort, SignalPort
        Edges: PowerBond, SignalBond, hasPowerPort, hasSignalPort
    TF : str
        The name of the Transformer

    Returns
    -------
    None

    side effects
    ------------
    The expression of the effort and flow of the PowerPorts and the PowerBond are updated if the energy is propagated through the Transformer
    
    """

    ports=getPowerPorts(G,TF)
    port_0=ports[0]
    port_1=ports[1]
    direction_0=checkBondDirection(G,port_0)
    direction_1=checkBondDirection(G,port_1)
    causality_0=G.nodes[port_0]['causality']
    causality_1=G.nodes[port_1]['causality']
    param_value = sympify(G.nodes[TF]['modelParameter']['param_1']['propertyValue'])
    if direction_0!=direction_1:
        if G.nodes[port_0][causality_0] is None and G.nodes[port_1][causality_0] is not None:
            G.nodes[port_0][causality_0]={}
            G.nodes[port_0][causality_0]['expression']=G.nodes[port_1][causality_0]['expression']*param_value
        if G.nodes[port_1][causality_1] is None and G.nodes[port_0][causality_1] is not None:
            G.nodes[port_1][causality_1]={}
            G.nodes[port_1][causality_1]['expression']=G.nodes[port_0][causality_1]['expression']*param_value
    else:
        if G.nodes[port_0][causality_0] is None and G.nodes[port_1][causality_0] is not None:
            G.nodes[port_0][causality_0]={}
            G.nodes[port_0][causality_0]['expression']=-G.nodes[port_1][causality_0]['expression']*param_value
        if G.nodes[port_1][causality_1] is None and G.nodes[port_0][causality_1] is not None:
            G.nodes[port_1][causality_1]={}
            G.nodes[port_1][causality_1]['expression']=-G.nodes[port_0][causality_1]['expression']*param_value
    
    powerBonds= [u for u in G.out_edges(port_0) if 'a' in G.edges[u].keys() and G.edges[u]['a']=='PowerBond']+ [u for u in G.in_edges(port_0) if 'a' in G.edges[u].keys() and G.edges[u]['a']=='PowerBond']
    viaPowerbond(G,powerBonds[0])
    powerBonds= [u for u in G.out_edges(port_1) if 'a' in G.edges[u].keys() and G.edges[u]['a']=='PowerBond']+ [u for u in G.in_edges(port_1) if 'a' in G.edges[u].keys() and G.edges[u]['a']=='PowerBond']
    viaPowerbond(G,powerBonds[0])
    return True

def viaGyrator(G,GY):
    """
    This function is used to propagate the energy through the Transformer

    Parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx
        Nodes: BondElement, JunctionStructure, PowerPort, SignalPort
        Edges: PowerBond, SignalBond, hasPowerPort, hasSignalPort
    TF : str
        The name of the Transformer

    Returns
    -------
    None

    side effects
    ------------
    The expression of the effort and flow of the PowerPorts and the PowerBond are updated if the energy is propagated through the Transformer
    
    """

    ports=getPowerPorts(G,GY)
    port_0=ports[0]
    port_1=ports[1]
    direction_0=checkBondDirection(G,port_0)
    direction_1=checkBondDirection(G,port_1)
    causality_0=G.nodes[port_0]['causality']
    causality_1=G.nodes[port_1]['causality']
    param_value = sympify(G.nodes[GY]['modelParameter']['param_1']['propertyValue'])
    if direction_0!=direction_1:
        if G.nodes[port_0][causality_0] is None and G.nodes[port_1][causality_1] is not None:
            G.nodes[port_0][causality_0]={}
            G.nodes[port_0][causality_0]['expression']=G.nodes[port_1][causality_1]['expression']*param_value
        if G.nodes[port_1][causality_1] is None and G.nodes[port_0][causality_0] is not None:
            G.nodes[port_1][causality_1]={}
            G.nodes[port_1][causality_1]['expression']=G.nodes[port_0][causality_0]['expression']*param_value
    else:
        if G.nodes[port_0][causality_0] is None and G.nodes[port_1][causality_1] is not None:
            G.nodes[port_0][causality_0]={}
            G.nodes[port_0][causality_0]['expression']=-G.nodes[port_1][causality_1]['expression']*param_value
        if G.nodes[port_1][causality_1] is None and G.nodes[port_0][causality_0] is not None:
            G.nodes[port_1][causality_1]={}
            G.nodes[port_1][causality_1]['expression']=-G.nodes[port_0][causality_0]['expression']*param_value
    
    powerBonds= [u for u in G.out_edges(port_0) if 'a' in G.edges[u].keys() and G.edges[u]['a']=='PowerBond']+ [u for u in G.in_edges(port_0) if 'a' in G.edges[u].keys() and G.edges[u]['a']=='PowerBond']
    viaPowerbond(G,powerBonds[0])
    powerBonds= [u for u in G.out_edges(port_1) if 'a' in G.edges[u].keys() and G.edges[u]['a']=='PowerBond']+ [u for u in G.in_edges(port_1) if 'a' in G.edges[u].keys() and G.edges[u]['a']=='PowerBond']
    viaPowerbond(G,powerBonds[0])
    return True

def checkBondDirection(G,port):
    if any(data.get('a') == 'PowerBond' for u, v, data in G.out_edges(port,data=True)):
        direction='out'
    elif any(data.get('a') == 'PowerBond' for u, v, data in G.in_edges(port,data=True)):
        direction='in'
    else:
        raise ValueError('The direction is not correct')
    return direction          

def nxBG_Energy(G):

    for edge in G.edges:
        if 'a' in G.edges[edge].keys() and G.edges[edge]['a']=='PowerBond':
            viaPowerbond(G,edge)
    for node in G.nodes:
        if G.nodes[node]['a']=='JunctionStructure' and G.nodes[node]['subClass']!='Transformer':
            viaJunction(G,node)
    for node in G.nodes:
        if G.nodes[node]['a']=='JunctionStructure' and G.nodes[node]['subClass']=='Transformer':
            viaTransformer(G,node)
    for node in G.nodes:
        if G.nodes[node]['a']=='JunctionStructure' and G.nodes[node]['subClass']!='Transformer':
            viaJunction(G,node)
    for node in G.nodes:
        if G.nodes[node]['a']=='JunctionStructure' and G.nodes[node]['subClass']=='Transformer':
            viaTransformer(G,node)
    # convert all the sympy expression to string
    for node in G.nodes:
        if G.nodes[node]['a']=='PowerPort':
            if G.nodes[node]['effort'] is not None:
                if 'expression' in G.nodes[node]['effort'].keys():
                    G.nodes[node]['effort']['expression']=str(G.nodes[node]['effort']['expression'])                    
            if G.nodes[node]['flow'] is not None:
                if 'expression' in G.nodes[node]['flow'].keys():
                    G.nodes[node]['flow']['expression']=str(G.nodes[node]['flow']['expression'])
    for edge in G.edges:
        if 'a' in G.edges[edge].keys() and G.edges[edge]['a']=='PowerBond':
            if G.edges[edge]['effort'] is not None:
                if 'expression' in G.edges[edge]['effort'].keys():
                    G.edges[edge]['effort']['expression']=str(G.edges[edge]['effort']['expression'])
            if G.edges[edge]['flow'] is not None:
                if 'expression' in G.edges[edge]['flow'].keys():
                    G.edges[edge]['flow']['expression']=str(G.edges[edge]['flow']['expression'])
    # check no None value in the effort and flow of the BondElement
    for node in G.nodes:
        if G.nodes[node]['a']=='PowerPort' and G.nodes[node]['isPortOf']=='BondElement':
            if G.nodes[node]['effort'] is None:
                raise ValueError(f'The effort of {node} is None')
            if G.nodes[node]['flow'] is None:
                raise ValueError(f'The flow of {node} is None')

def nxBG_preCellML(G):
    input_vars=[]
    output_vars=[]
    conservation_laws={}
    for node in G.nodes:
        if G.nodes[node]['a']=='PowerPort' and G.nodes[node]['isPortOf']=='BondElement':
            if G.nodes[node]['causality']=='flow':
                input_vars.append(G.nodes[node]['flow']['propertyName'])
                output_var='e_'+node
                output_vars.append(output_var)
                conservation_laws.update({output_var:G.nodes[node]['effort']['expression']})
            elif G.nodes[node]['causality']=='effort':
                input_vars.append(G.nodes[node]['effort']['propertyName'])
                output_var='f_'+node
                output_vars.append(output_var)
                conservation_laws.update({output_var:G.nodes[node]['flow']['expression']})
    return input_vars, output_vars, conservation_laws        

def save_nxBG_json(G, filename='nx_BG.json'):
    with open(filename, 'w') as f:
        json.dump(json_graph.node_link_data(G,edges="edges"), f, indent=4)

def load_nxBG_json(filename):
    with open(filename, 'r') as f:
        data = json.load(f)
    G = json_graph.node_link_graph(data, edges="edges", directed=True)
    return G

def pyvis_BG(G):
    pyvis_net = Network(cdn_resources='in_line', directed=True)
    pyvis_net.from_nx(G)
    return pyvis_net

def pyvis_net_fomat(pyvis_net):
    """
    This function is used to format the bond graph (pyvis_net) for visualization
        
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
                elif node['subClass']=='Transformer':
                    node['shape']='square'
                    node['color']={'background':'yellow'}
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



def nxBG_pyvis_net(pyvis_net, filename='nx_BG.html'):
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

def save_nxBG_pyvis_net(G, filename='nx_BG.html'):

    pyvis_net = pyvis_BG(G)
    pyvis_net_fomat(pyvis_net)
    nxBG_pyvis_net(pyvis_net, filename)


if __name__ == "__main__": 
    import webbrowser
    import os
    G=nx.DiGraph()
    # Get the default web browser
    browser = webbrowser.get()
    path_='./data/'
    fmatrix='./data/SLC5_f.csv'
    rmatrix='./data/SLC5_r.csv'
    nxBG(G, fmatrix, direction='e2f', connection='PowerBond')
    nxBG(G, rmatrix, direction='f2e', connection='PowerBond')    
    
    save_nxBG_json(G, path_+'nx_BG_0.json')
    save_nxBG_pyvis_net(G, path_+'nx_BG_0.html')
    browser.open(os.path.abspath(path_ + 'nx_BG_0.html'))   
    
    nxBG_addJunctions(G)

    nxBG_initEnergy(G)

    nxBG_Energy(G)

    save_nxBG_json(G, path_+'nx_BG_Energy.json')
    save_nxBG_pyvis_net(G, path_+'nx_BG_Energy.html')
    
    browser.open(os.path.abspath(path_ + 'nx_BG_Energy.html')) 
    G_1=load_nxBG_json(path_+'nx_BG_Energy_1.json') 

    G2=nx.compose(G,G_1)
    nxBG_addJunction(G2)
    save_nxBG_json(G2, path_+'nx_BG_Energy_2.json')
    save_nxBG_pyvis_net(G2, path_+'nx_BG_Energy_2.html')
    browser.open(os.path.abspath(path_ + 'nx_BG_Energy_2.html'))

    input_vars, output_vars, conservation_laws=nxBG_preCellML(G)
    print('Input variables:', input_vars)
    print('Output variables:', output_vars)
    print('Conservation laws:')
    for output_var, expression in conservation_laws.items():
        print(f'{output_var} = {expression}')
    