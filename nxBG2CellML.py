def nxBG_refine_subclass(G):
    """
    This function is used to specify the subclass of the BondElement based on the user input
    
    """
    BondElement_subclasses = ['Dissipative', 'Storage', 'Source', 'Unknown']
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
                    
def _nxBG_orientation_out(G):
    """
    This function is used to specify the orientation of the power ports of the BondElement based on the user input

    Parameters
    ----------
    G : nx.DiGraph
        The bond graph built using networkx
        Nodes: BondElement, JunctionStructure, PowerPort, SignalPort
        Edges: PowerBond, SignalBond, hasPowerPort, hasSignalPort

    Returns
    -------
    None

    side effects
    ------------
    The orientation of the power ports of the BondElement are updated based on the user input
    
    """
    PowerPorts=[node for node in G.nodes if G.nodes[node]['a']=='PowerPort' and G.nodes[node]['isPortOf']=='BondElement']
    for idx, node in enumerate(PowerPorts):
        print(f"{idx}. {node}")
    print("Select the power ports whose orientation should be 'out' (the default is 'in'):")
    print("Enter the number of the power ports separated by comma:")
    out_ports_num = input().split(',')
    for port_num in out_ports_num:
        G.nodes[PowerPorts[int(port_num)]]['orientation'] = 'out'

def nxBG_component(G,BG_component):
    bondElements=[node for node in G.nodes if G.nodes[node]['a']=='BondElement']
    
    for node in bondElements:
        if G.nodes[node]['component'] in BG_component.keys():
            powerPorts=_getPowerPorts(G,node)
            for port in powerPorts:
                port_num=port.split('_')[-1]
                G.nodes[port]['orientation'] = BG_component[G.nodes[node]['component']]['ports'][port_num]['orientation']
                G.nodes[port]['effort']['propertyName'] =BG_component[G.nodes[node]['component']]['ports'][port_num]['effort']['propertyName']
                G.nodes[port]['flow']['propertyName'] = BG_component[G.nodes[node]['component']]['ports'][port_num]['flow']['propertyName']
                G.nodes[port]['effort']['units'] = BG_component[G.nodes[node]['component']]['ports'][port_num]['effort']['units']
                G.nodes[port]['flow']['units'] = BG_component[G.nodes[node]['component']]['ports'][port_num]['flow']['units']

            signalPorts=_getSignalPorts(G,node)
            
            for port in signalPorts:
                port_num=port.split('_')[-1]
                G.nodes[port]['signal']['propertyName'] = BG_component[G.nodes[node]['component']]['ports'][port_num]['signal']['propertyName']
                G.nodes[port]['signal']['units'] = BG_component[G.nodes[node]['component']]['ports'][port_num]['signal']['units']

            if 'params' in BG_component[G.nodes[node]['component']].keys():
                G.nodes[node]['modelParameter'] = {}
                for param in BG_component[G.nodes[node]['component']]['params']:
                    G.nodes[node]['modelParameter'][param] = {}
                    G.nodes[node]['modelParameter'][param]['propertyName'] = BG_component[G.nodes[node]['component']]['params'][param]['propertyName']
                    G.nodes[node]['modelParameter'][param]['units'] = BG_component[G.nodes[node]['component']]['params'][param]['units']
                    G.nodes[node]['modelParameter'][param]['propertyValue'] = BG_component[G.nodes[node]['component']]['params'][param]['propertyValue']
            
            if 'state_vars' in BG_component[G.nodes[node]['component']].keys():
                G.nodes[node]['modelState'] = {}
                for state_var in BG_component[G.nodes[node]['component']]['state_vars']:
                    G.nodes[node]['modelState'][state_var] = {}
                    G.nodes[node]['modelState'][state_var]['propertyName'] = BG_component[G.nodes[node]['component']]['state_vars'][state_var]['propertyName']
                    G.nodes[node]['modelState'][state_var]['units'] = BG_component[G.nodes[node]['component']]['state_vars'][state_var]['units']
                    G.nodes[node]['modelState'][state_var]['propertyValue'] = BG_component[G.nodes[node]['component']]['state_vars'][state_var]['propertyValue']
            
            if 'constitutive_relations' in BG_component[G.nodes[node]['component']].keys():
                G.nodes[node]['hasConstitutiveRelation'] = BG_component[G.nodes[node]['component']]['constitutive_relations']