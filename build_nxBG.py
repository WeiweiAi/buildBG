import networkx as nx
from buildBG import load_matrix
import matplotlib.pyplot as plt
from networkx.readwrite import json_graph
import json
from pyvis.network import Network
from sympy import *

def nxBG(G, matrix, direction='e2f', connection='PowerBond'):
    """
    This function is used to build the bond graph using networkx based on stoichiometry-like matrix
    Note: The 0 or 1 junctions may not be added to the bond graph, which can be added by nxBG_addJunction function
    
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
    cName=eName+fName
    cPort=ePort+fPort
    cID=eID+fID
    for i in range(len(cName)):
        if cName[i] not in G.nodes:
            if cID[i]=='1': # BondGraph Junction
                G.add_node(cName[i], a='JunctionStructure', subClass='OneJunctionStructure')
            elif cID[i]=='0': # BondGraph Junction
                G.add_node(cName[i], a='JunctionStructure', subClass='ZeroJunctionStructure')
            else: # BondGraph Element
                G.add_node(cName[i], a='BondElement')
        iport=f'{cName[i]}_{cPort[i]}'
        if connection=='PowerBond':
            if i<len(eName):
                if cID[i]=='in':
                    G.add_node(iport, a='PowerPort', isPortOf=G.nodes[cName[i]]['a'], effort=None, flow=None, causality='effort-out', positive='in')
                elif cID[i]=='out':
                    G.add_node(iport, a='PowerPort', isPortOf=G.nodes[cName[i]]['a'], effort=None, flow=None, causality='effort-out', positive='out')
                else:
                    G.add_node(iport, a='PowerPort', isPortOf=G.nodes[cName[i]]['a'], effort=None, flow=None, causality='effort-out', positive='unknown')
            if i>=len(eName):
                if cID[i]=='in':
                    G.add_node(iport, a='PowerPort', isPortOf=G.nodes[cName[i]]['a'], effort=None, flow=None, causality='flow-out', positive='in')
                elif cID[i]=='out':
                    G.add_node(iport, a='PowerPort', isPortOf=G.nodes[cName[i]]['a'], effort=None, flow=None, causality='flow-out', positive='out')
                else:
                    G.add_node(iport, a='PowerPort', isPortOf=G.nodes[cName[i]]['a'], effort=None, flow=None, causality='flow-out', positive='unknown')
            
            G.add_edge(cName[i],iport, relationship='hasPowerPort')
        elif connection=='SignalBond':
            G.add_node(iport, a='SignalPort', isPortOf=G.nodes[cName[i]]['a'], signal=None)
            G.add_edge(cName[i],iport, relationship='hasSignalPort')
        else:
            raise ValueError('The connection type is not correct')
    for i in range(len(eName)):
        for j in range(len(fName)):
            iport=f'{eName[i]}_{ePort[i]}'
            jport=f'{fName[j]}_{fPort[j]}'
            if N[i][j]>0:
                G.add_edge(iport,jport, weight=str(N[i][j]), a=connection, effort=None, flow=None)
            elif N[i][j]<0:   
                G.add_edge(jport,iport, weight=str(-N[i][j]), a=connection, effort=None, flow=None)
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
                G.add_edge(jport,edge[1], weight=G.edges[edge]['weight'], a=G.edges[edge]['a'], effort=G.edges[edge]['effort'], flow=G.edges[edge]['flow'])
                G.remove_edge(edge[0],edge[1])
            for edge in in_edges:
                jPort+=1
                jport=f'{jName}_{jPort}'
                G.add_node(jport, a='PowerPort', isPortOf=G.nodes[jName]['a'], effort=None, flow=None, causality=G.nodes[node]['causality'])
                G.add_edge(jName,jport, relationship='hasPowerPort')
                G.add_edge(edge[0],jport, weight=G.edges[edge]['weight'], a=G.edges[edge]['a'], effort=G.edges[edge]['effort'], flow=G.edges[edge]['flow'])
                G.remove_edge(edge[0],edge[1]) 
            jPort=0
            jport=f'{jName}_{jPort}'
            if G.nodes[node]['causality']=='effort-out':
                G.nodes[jName]['subClass']='ZeroJunctionStructure'
                G.add_node(jport, a='PowerPort', isPortOf=G.nodes[jName]['a'], effort=None, flow=None, causality='flow-out')
                if G.nodes[node]['positive']=='in':
                    G.add_edge(jport,node, weight='1', a='PowerBond', effort=None, flow=None)
                elif G.nodes[node]['positive']=='out':
                    G.add_edge(node,jport, weight='1', a='PowerBond', effort=None, flow=None)
                else:
                    raise ValueError('The positive flow direction is unknown')
                
            elif G.nodes[node]['causality']=='flow-out':
                G.nodes[jName]['subClass']='OneJunctionStructure'
                G.add_node(jport, a='PowerPort', isPortOf=G.nodes[jName]['a'], effort=None, flow=None, causality='effort-out')
                if G.nodes[node]['positive']=='in':
                    G.add_edge(jport,node,weight='1', a='PowerBond', effort=None, flow=None)
                elif G.nodes[node]['positive']=='out':
                    G.add_edge(node,jport, weight='1', a='PowerBond', effort=None, flow=None)
                else:
                    raise ValueError('The positive flow direction is unknown')
            else:
                raise ValueError('The causality is not correct')      
            G.add_edge(jName,jport, relationship='hasPowerPort')

def nxBG_addJunction_multi(G):
    """
    This function is used to replace the PowerBond with the Transformer in case the PowerBond has the weight not equal to 1
    Note: The Transformer is a subclass of the BondElement, which is different from https://celldl.org/ontologies/bondgraph#MultiplierJunctionStructure

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
    PowerBonds=[edge for edge in G.edges if 'a' in G.edges[edge].keys() and  G.edges[edge]['a']=='PowerBond' and float(G.edges[edge]['weight'])!=1]    
    for edge in PowerBonds:
        source=edge[0]
        target=edge[1]
        TFName='TF_'+source+'_'+target
        G.add_node(TFName, a='BondElement', subClass='Transformer', modelParameter={'param_1':{'value':float(G.edges[edge]['weight'])}})
        # add the power ports of the transformer
        port_0=f'{TFName}_0'
        port_1=f'{TFName}_1'
        G.add_node(port_0, a='PowerPort', isPortOf=G.nodes[TFName]['a'], effort=None, flow=None, causality='flow-out', positive='in')
        G.add_edge(TFName,port_0, relationship='hasPowerPort')
        G.add_node(port_1, a='PowerPort', isPortOf=G.nodes[TFName]['a'], effort=None, flow=None, causality='effort-out', positive='out')
        G.add_edge(TFName,port_1, relationship='hasPowerPort')
        if G.nodes[source]['causality']=='effort-out':
            G.add_edge(source,port_0, weight='1', a='PowerBond', effort=None, flow=None)
            G.add_edge(port_1,target, weight='1', a='PowerBond', effort=None, flow=None)
        elif G.nodes[source]['causality']=='flow-out':
            G.add_edge(source,port_1, weight='1', a='PowerBond',  effort=None, flow=None)
            G.add_edge(port_0,target, weight='1', a='PowerBond',  effort=None, flow=None)
        else:
            raise ValueError('The causality is not correct')
        G.remove_edge(source,target)

def nxBG_initEnergy(G):
    """
    This function is used to initialize:
    the output variable of the PowerPort of the BondElement based on the causality of the PowerBond, either effort or flow

    """
    for node in G.nodes:
        if G.nodes[node]['a']=='PowerPort' and G.nodes[node]['isPortOf']=='BondElement':
            if G.nodes[node]['causality']=='effort-out' and G.nodes[node]['effort'] is None:
                G.nodes[node]['effort']='e_'+ node
            if G.nodes[node]['causality']=='flow-out' and G.nodes[node]['flow'] is None:
                G.nodes[node]['flow']='f_'+ node

def PowerPort2PowerBond(G,port):
    """
    This function is used to propagate the effort and flow from the power port to the power bond    
    """
    power_bonds=[]
    if G.nodes[port]['a']=='PowerPort':
        power_bonds+=[(u, v) for u, v, data in G.out_edges(port,data=True) if data.get('a') == 'PowerBond']
        power_bonds+=[(u, v) for u, v, data in G.in_edges(port,data=True) if data.get('a') == 'PowerBond']
        for power_bond in power_bonds:
            if G.nodes[port]['effort'] is not None:
                G.edges[power_bond]['effort']=G.nodes[port]['effort']
            if G.nodes[port]['flow'] is not None:
                G.edges[power_bond]['flow']=G.nodes[port]['flow']

def PowerBond2PowerPort(G,edge):
    """
    This function is used to propagate the effort and flow from the power bond to the power port
    """
    if 'a' in G.edges[edge].keys() and G.edges[edge]['a'] == 'PowerBond':
       source=edge[0]
       target=edge[1]
       if G.nodes[source]['effort'] is None and G.edges[edge]['effort'] is not None:
            G.nodes[source]['effort']=G.edges[edge]['effort']
       if G.nodes[source]['flow'] is None and G.edges[edge]['flow'] is not None:
            G.nodes[source]['flow']=G.edges[edge]['flow']
       if G.nodes[target]['effort'] is None and G.edges[edge]['effort'] is not None:
            G.nodes[target]['effort']=G.edges[edge]['effort']
       if G.nodes[target]['flow'] is None and G.edges[edge]['flow'] is not None:
            G.nodes[target]['flow']=G.edges[edge]['flow']

def _checkJunction_sharedVar(G,nodeJunction):
    """
    This function is used to check the shared variable of the junction

    parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx
        Nodes: BondElement, JunctionStructure, PowerPort, SignalPort
        Edges: PowerBond, SignalBond, hasPowerPort, hasSignalPort
    nodeJunction : str
        The name of the junction node

    Returns
    -------
    shared : str
        The shared variable of the junction
    shared_var : str
        The shared variable type, either effort or flow
    solved : bool
        The flag to indicate if the shared variable is solved
        i.e., all the power ports of the junction have the same shared variable
    ports : list
        The list of the power ports of the junction

    """
    if G.nodes[nodeJunction]['subClass']=='ZeroJunctionStructure':
        shared_var='effort'
    elif G.nodes[nodeJunction]['subClass']=='OneJunctionStructure':
        shared_var='flow'
    ports=[out_edge[1] for out_edge in G.out_edges(nodeJunction) if 'relationship' in G.edges[out_edge]] # get the power ports of the junction
    shared=None
    nports_known=0
    for port in ports:
        if G.nodes[port][shared_var] is not None:
            shared=G.nodes[port][shared_var] # get the shared variable
            nports_known+=1
    if nports_known==len(ports):
        PowerPort2PowerBond(G,port) # may not be necessary
        solved=True
    else:
        solved=False
    return shared, shared_var, solved,ports

def _updateJunction_sharedVar(G,nodeJunction):
    """
    This function is used to:
      update the shared variable of the ports of the junction;
      update the shared variable of the power bond connected to the junction;
      update the shared variable of the power port of the other side of the power bond;
      propagate the shared variable to the other junctions connected to the junction

      Note: This function is critical and needs to be tested carefully
    """

    shared, shared_var, solved,ports=_checkJunction_sharedVar(G,nodeJunction)
    if (shared is not None) and (not solved):
        for port in ports:
            G.nodes[port][shared_var]=shared # update the shared variable for the power ports
            edges=[(u, v) for u, v, data in G.out_edges(port,data=True) if data.get('a') == 'PowerBond']
            if len(edges)==1:
                port_v=edges[0][1]
            else:
                edges+=[(u, v) for u, v, data in G.in_edges(port,data=True) if data.get('a') == 'PowerBond']
                if len(edges)==1:
                    port_v=edges[0][0]
            if len(edges)==1:
                G.edges[edges[0]][shared_var]=shared # update the shared variable for the power bond
                G.nodes[port_v][shared_var]=shared # update the shared variable for the power port of the other side of the power bond
                for u, v, data in G.in_edges(port_v,data=True) :
                    if data.get('relationship') == 'hasPowerPort':
                        node_v=u
                        if G.nodes[node_v]['a']=='JunctionStructure':
                            _updateJunction_sharedVar(G,node_v)
                            break              
            else:
                raise ValueError('The number of edges is not correct')
        return True
    else:
        return False

def updateJunction_sharedVar(G):        
    junctions=[node for node in G.nodes if G.nodes[node]['a']=='JunctionStructure']
    for nodeJunction in junctions:
        _updateJunction_sharedVar(G,nodeJunction)

def checkBondDirection(G,port):
    if any(data.get('a') == 'PowerBond' for u, v, data in G.out_edges(port,data=True)):
        direction='out'
    elif any(data.get('a') == 'PowerBond' for u, v, data in G.in_edges(port,data=True)):
        direction='in'
    else:
        raise ValueError('The direction is not correct')
    return direction

def updateJunction_conservation(G):
    """
    This function is used to:
      update the effort or flow of the power port of the junction based on the conservation law
      propagate the effort or flow to the power bond connected to the junction

    """

    junctions=[node for node in G.nodes if G.nodes[node]['a']=='JunctionStructure']
    for nodeJunction in junctions:
        shared, shared_var, solved,ports=_checkJunction_sharedVar(G,nodeJunction)
        if shared_var=='effort':
            compute_var='flow'
        elif shared_var=='flow':
            compute_var='effort'
        else:
            raise ValueError('The shared variable is not correct')
        compute_var_init=0
        for port in ports:
            if G.nodes[port][compute_var] is None:
                compute_port=port
                compute_direction=checkBondDirection(G,port)
                break
        for port in ports:
            if G.nodes[port][compute_var] is not None:
                if checkBondDirection(G,port)==compute_direction:
                    compute_var_init-=symbols(G.nodes[port][compute_var])
                else:
                    compute_var_init+=symbols(G.nodes[port][compute_var])
        G.nodes[compute_port][compute_var]=ccode(compute_var_init)
        PowerPort2PowerBond(G,compute_port)
        edges=[(u, v) for u, v, data in G.out_edges(compute_port,data=True) if data.get('a') == 'PowerBond']
        edges+=[(u, v) for u, v, data in G.in_edges(compute_port,data=True) if data.get('a') == 'PowerBond']
        for edge in edges:
            PowerBond2PowerPort(G,edge)           

def nxBG_refine_subclass(G):
    """
    This function is used to specify the subclass of the BondElement based on the user input
    
    """
    BondElement_subclasses = ['Dissipative', 'Storage', 'Source', 'Unknown','Transformer','Gyrator']
    # print the subclass with the numberings
    for idx, subclass in enumerate(BondElement_subclasses, start=1):
        print(f"{idx}. {subclass}")
    for node in G.nodes:
        if G.nodes[node]['a']=='BondElement' and 'subClass' not in G.nodes[node].keys():
            print(f"Refining node: {node}")
            print("Select a subclass:")
            subclass = input()
            G.nodes[node]['subClass'] = BondElement_subclasses[int(subclass)-1]

def prelink_nxBG_CellML(G):
    """
    This function is used to specify the number of modelParameter and stateVariable and add them to the BondElement based on the user input
    """
    for node in G.nodes:
        if G.nodes[node]['a'] == 'BondElement':
            # ask users if there is any parameter for this node
            print('How many parameters for this node?')
            param = input()
            if int(param)>0:
                G.nodes[node]['modelParameter'] = {}
                for p in range(int(param)):
                    G.nodes[node]['modelParameter'][f'param_{p+1}'] = None               
            # ask users if there is any state variable for this node
            print('How many state variables for this node?')
            state = input()
            if int(state)>0:
                G.nodes[node]['stateVariable'] = {}
                for s in range(int(state)):
                    G.nodes[node]['stateVariable'][f'q_{s+1}']= None

def nxBG_addJunctions(G):
    """
    This function is used to add transformers and 0 or 1 junctions
    """
    nxBG_addJunction_multi(G)
    nxBG_addJunction(G)


def nxBG_propogateEnergy(G):

    for node in G.nodes:
        PowerPort2PowerBond(G,node)

    for edge in G.edges:
        PowerBond2PowerPort(G,edge)

    updateJunction_sharedVar(G)

    updateJunction_conservation(G)
 
    updateJunction_sharedVar(G)

    for edge in G.edges:
        PowerBond2PowerPort(G,edge)

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
    conservation_laws=[]
    for node in G.nodes:
        if G.nodes[node]['a']=='PowerPort' and G.nodes[node]['isPortOf']=='BondElement':
            if G.nodes[node]['causality']=='flow-out':
                input_vars.append(G.nodes[node]['flow'])
                output_var='e_'+node
                output_vars.append(output_var)
                conservation_laws.append(output_var+'='+G.nodes[node]['effort'])
            elif G.nodes[node]['causality']=='effort-out':
                input_vars.append(G.nodes[node]['effort'])
                output_var='f_'+node
                output_vars.append(output_var)
                conservation_laws.append(output_var+'='+G.nodes[node]['flow'])
    return input_vars, output_vars, conservation_laws        

if __name__ == "__main__": 

    G=nx.DiGraph()
    # Create a PyVis network with cdn_resources set to 'in_line'
    pyvis_net = Network(cdn_resources='in_line',directed =True)
    fmatrix='./data/SLC5_f.csv'
    rmatrix='./data/SLC5_r.csv'
    nxBG(G, fmatrix, direction='e2f', connection='PowerBond')
    nxBG(G, rmatrix, direction='f2e', connection='PowerBond')    
    nxBG_addJunctions(G)    
    nxBG_initEnergy(G)
    nxBG_propogateEnergy(G)
    with open('nx_BG_propogateEnergy.json', 'w') as f:
        json.dump(json_graph.node_link_data(G,edges="edges"), f, indent=4)
    
    pyvis_net.from_nx(G)

    # Save the network to an HTML file with utf-8 encoding
    html_file = "BG.html"
    with open(html_file, 'w', encoding='utf-8') as f:
        f.write(pyvis_net.generate_html())

    # Display the network in the default web browser
    import webbrowser
    webbrowser.open(html_file)

    input_vars, output_vars, conservation_laws=nxBG_preCellML(G)
    print('Input variables:', input_vars)
    print('Output variables:', output_vars)
    print('Conservation laws:')
    for law in conservation_laws:
        print(law)
