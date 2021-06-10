import networkx as nx
import numpy as np
from scipy.cluster import hierarchy
from ChromProcess import build_cluster_tree

def graph_from_linkage(linkage_mat, id_modifier = ''):
    '''
    linkage_mat: scipy condensed linkage matrix
        linkage matrix

    id_modifier: str
        name modifier (can be used to distinguish different clusterings)
    '''
    rootnode, nodelist = hierarchy.to_tree(linkage_mat, rd= True)

    G = nx.DiGraph()

    for n in nodelist:
        source = '{}{}'.format(n.id, id_modifier)
        G.add_node(source, distance = n.dist, leaf = n.is_leaf(),
                    id = n.id)
        if n.count > 1:
            target_left = '{}{}'.format(n.left.id, id_modifier)

            G.add_node(target_left, distance = n.left.dist, leaf = n.left.is_leaf(),
                        id = n.left.id)
            G.add_edge(source, target_left, weight = abs(n.left.dist - n.dist),
                        dist = abs(n.left.dist - n.dist))

            target_right = '{}{}'.format(n.right.id, id_modifier)
            G.add_node(target_right, distance = n.right.dist, leaf = n.right.is_leaf(),
                        id = n.right.id)
            G.add_edge(source, target_right, weight = abs(n.right.dist - n.dist),
                        dist = abs(n.right.dist - n.dist))

    return G

def graphviz_layout_networkx(network, render_engine = 'fdp'):
    '''
    Uses graphviz to generate a layout from a networkx graph.

    Layouts from the graphviz documentation:
    dot − filter for drawing directed graphs
    neato − filter for drawing undirected graphs
    twopi − filter for radial layouts of graphs
    circo − filter for circular layout of graphs
    fdp − filter for drawing undirected graphs
    sfdp − filter for drawing large undirected graphs
    patchwork − filter for squarified tree maps
    osage − filter for array-based layouts

    Parameters
    ----------
    network: Networkx DiGraph
        Network to be represented.
    render_engine: str
        Layout render engine for graphviz.
    Returns
    -------
    pos: dict
        {node:[float(x),float(y)]}
    '''
    import json
    from graphviz import Digraph

    # Create a graph with graphviz to plot a scheme of the network
    dot = Digraph(comment = '',
                  engine = render_engine,
                  strict = 'True',
                  format = 'json')

    for n in network.nodes: # add in nodes
        dot.node(n,n)

    for e in network.edges: # Create edges between the nodes.
        dot.edge(e[0],e[1])

    json_string = dot.pipe().decode()

    y = json.loads(json_string)

    pos = {}
    for o in y['objects']:
        pos[o['name']] = [float(x) for x in o['pos'].split(',')]

    return pos

def set_network_coords(network, pos):
    import numpy as np

    xmin = np.mean([pos[p][0] for p in pos])
    ymin = np.mean([pos[p][1] for p in pos])
    for n in network.nodes:
        if n in pos:
            network.nodes[n]['pos'] = pos[n]
        else:
            network.nodes[n]['pos'] = (xmin, ymin)

    return network

def get_network_lineplot(G):
    '''
    Parameters
    ----------
    G: networkx DiGraph
        Graph to extract nodes from.

    Returns
    -------
    net_lines: numpy 2D array
        Coordinates for plotting a line plot of the network.
    '''

    import numpy as np

    net_lines = []
    for e in G.edges:
        for n in e:
            net_lines.append(G.nodes[n]["pos"])
        net_lines.append((np.nan,np.nan))

    net_lines = np.array(net_lines)
    net_lines = net_lines.T
    return net_lines

def get_network_scatter(G):
    '''
    Parameters
    ----------
    G: networkx DiGraph
        Graph to extract nodes from.

    Returns
    -------
    net_lines: numpy 2D array
        Coordinates for plotting a line plot of the network.
    '''

    import numpy as np

    net_scatter = []
    for n in G.nodes:
        net_scatter.append(G.nodes[n]["pos"])

    net_scatter = np.array(net_scatter)
    net_scatter = net_scatter.T
    return net_scatter

def normalise_network_coordinates(G):
    '''
    Parameters
    ----------
    G: networkx DiGraph
        Graph to extract nodes from.

    Returns
    -------
    None
    '''
    import numpy as np
    coords = get_network_scatter(G)

    net_width = (np.max(coords[0])-np.min(coords[0]))
    net_height = (np.max(coords[1])-np.min(coords[1]))

    if net_width == 0:
        net_width = np.max(coords[0])
    if net_height == 0:
        net_height = np.max(coords[1])
        
    for n in G.nodes:
        pos = G.nodes[n]['pos']
        a = pos[0]/net_width
        b = pos[1]/net_height
        G.nodes[n]['pos'] = (a,b)
