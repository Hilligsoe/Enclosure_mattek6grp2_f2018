# -*- coding: utf-8 -*-
"""
Created on Thu Sep 07 14:46:11 2017

@author: cht15
"""
import time
import numpy as np
from itertools import combinations
import networkx as nx
from numba import jit


def _slow_edges(G, radius, p=2):
    """Returns edge list of node pairs within `radius` of each other
       using Minkowski distance metric `p`

    Works without scipy, but in `O(n^2)` time.
    """
    # TODO This can be parallelized.
    # TODO This we might be able to do it with the fast edges - 'O(n)'
    edges = []
    d = {}
    d4 = {}
    for (u, pu), (v, pv) in combinations(G.nodes(data='pos'), 2):
        distance = sum(abs(a - b) ** p for a, b in zip(pu, pv))
        if distance <= radius ** p:
            edges.append((u, v))
            d[(u, v)] = np.sqrt(distance)
            d4[(u, v)] = distance**2
    return edges, d, d4


class GraphCreation2D:
    def __init__(self, x1, x2, y1, y2, intensity):
        # init Height, width and intensity
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        self.intensity = intensity

    def __call__(self, Tr, SD={'S': (-5, 0), 'D': (5, 0)},
                 Tv={'R': 2, 'packagesize': 1, 'Processinglag': 1}):

        N = np.random.poisson(self.intensity*(np.abs(self.x1 - self.x2)) *
                              (np.abs(self.y1-self.y2)))
        # Number of points by poisson
#        t = time.time()
        location = {i: (np.random.uniform(self.x1, self.x2),
                        np.random.uniform(self.y1, self.y2)) for i in range(N)}
#        print('Time to create %.f nodes: ' % N + repr(time.time()-t))
        Graph = nx.Graph()
        Graph.add_nodes_from(range(N))
        Graph.add_nodes_from(SD)
        for key in SD:
            Graph.node[key]['pos'] = SD[key]
            location[key] = SD[key]
        nx.set_node_attributes(Graph, location, 'pos')
#        t = time.time()
        edges, distance, d4 = _slow_edges(Graph, Tr)
#        print('Time to create edges: ' + repr(time.time()-t) + " seconds")
        Graph.add_edges_from(edges)
        nx.set_edge_attributes(Graph, distance, 'distance')
        nx.set_edge_attributes(Graph, d4, 'd4')
        nx.set_edge_attributes(Graph, 0, 'hop_limit')
        return Graph, location, N, distance


@jit
def _gaincalculation(width, efficiency):
    g1 = 2/(1 - np.cos(width/(efficiency*2.)))
    Delta = 2/(1 - np.cos(width/2))
    g2 = (Delta - g1)/(Delta - 1)
    return g1, g2


@jit
def sufficient_range(intensity, f):
    return np.sqrt((np.log(intensity)+f(intensity))/(np.pi*intensity))


@jit
def diskrange(Systemvariable={'width': np.radians(30), 'efficiency': 0.8,
                              'R': 2, 'power': 1, 'noise': 1}):
    #  width, efficiency, R, power, noise
    G1, g2 = _gaincalculation(Systemvariable['width'],
                              Systemvariable['efficiency'])
    peak = ((2**Systemvariable['R']) - 1.)*Systemvariable['noise']
    return (G1*Systemvariable['power']/(peak))**(0.25), G1, g2


if __name__ == '__main__':
    np.random.seed(3)

    systemvariables = {'width': np.radians(30), 'efficiency': 0.8, 'R': 2,
                       'power': 1, 'noise': 1}

    Transmissionvariables = {'R': systemvariables['R'], 'packagesize': 1,
                             'Processinglag': 1}

    Transmissionrange, P_main, P_side = diskrange(systemvariables)

    ex = GraphCreation2D(x1=-10, x2=10, y1=-10, y2=10, intensity=4)

    Graph, Location, N, D = ex(Tr=Transmissionrange, Tv=Transmissionvariables)
    colour = np.array(['C0' for i in range(len(Graph)-2)])
    colour = np.hstack((colour, ['C1', 'C1']))
    nx.draw_networkx(Graph, Location, with_labels=False, node_size=20,
                     node_color=colour, edge_color='slategrey')
    print('Is the graph fully connected?: ' + repr(nx.is_connected(Graph)))
    if nx.has_path(Graph, 'S', 'D'):
        path0 = nx.shortest_path(Graph, source='S', target='D')
        path1 = nx.shortest_path(Graph, source='S', target='D',
                                 weight='distance')
        path2 = nx.shortest_path(Graph, source='S', target='D',
                                 weight='d4')
        path3 = nx.shortest_path(Graph, source='S', target='D',
                                 weight=('d4', 'hop_limit'))

        print('The shortest path between Source and Destination is \n'
              + repr(path0) + '\n weighted with distance \n:' + repr(path1))

        route_edges0 = [(path0[n], path0[n+1]) for n in range(len(path0)-1)]
        route_edges1 = [(path1[n], path1[n+1]) for n in range(len(path1)-1)]
        route_edges2 = [(path2[n], path2[n+1]) for n in range(len(path2)-1)]
        route_edges3 = [(path3[n], path3[n+1]) for n in range(len(path3)-1)]
        nx.draw_networkx_edges(Graph, Location, route_edges0, edge_color='C2',
                               width=4)
        nx.draw_networkx_edges(Graph, Location, route_edges1, edge_color='C3',
                               width=4)
        nx.draw_networkx_edges(Graph, Location, route_edges2, edge_color='C4',
                               width=4)
        nx.draw_networkx_edges(Graph, Location, route_edges3, edge_color='C5',
                               width=4)
    else:
        print("There isn't a path between Source and Destination")

    print('\n With a Transmissions range of ' + repr(Transmissionrange) +
          ' and an Transmissionrate of ' + repr(systemvariables['R']))
