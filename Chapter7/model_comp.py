# -*- coding: utf-8 -*-
"""
Created on Thu May 03 15:56:31 2018

@author: MATTEK6 grp 2
"""

import numpy as np
import networkx as nx
from path_import import diskrange, GraphCreation2D, _gaincalculation
from wmn_main import system_transmission, Rate, _interferencelimit
import scipy as sc
import warnings
warnings.filterwarnings("ignore")


if __name__ == '__main__':
    # initializing varibles for summation of results to determine averge
    abs0 = 0
    pro0 = 0
    abs1 = 0
    pro1 = 0
    l_abs0 = 0
    l_pro0 = 0
    l_abs1 = 0
    l_pro1 = 0
    throughput_sum0 = 0
    throughput_sum1 = 0
    latency_sum0 = 0
    latency_sum1 = 0
    ET0_sum = 0
    ET1_sum = 0
    EL0_sum = 0
    EL1_sum = 0
    connection_sum = 0
    above0 = 0
    above1 = 0

    # set random seed (3 is used for documented results)
    np.random.seed(3)

    # intital system settings
    intensity = 5

    width = np.radians(30)
    efficiency = 0.8
    P_main, P_side = _gaincalculation(width, efficiency)
    print('Gain: %.5f, %.5f' % (P_main, P_side))
    systemvariables = {'width': width, 'efficiency': efficiency,
                       'R': Rate(intensity, P_main), 'power': 1, 'noise': 1}
    # SourceDestination = {'S': (0, 0), 'D': (1, 0)}
    # SourceDestination = {'S': (-0.25, 0), 'D': (.25, 0)}
    SourceDestination = {'S': (-2.5, 0), 'D': (2.5, 0)}
    Transmissionvariables = {'R': systemvariables['R'], 'packagesize': 1,
                             'Processinglag': 1, 'intensity': intensity}
    p_trans = 0.5
    maxiter = 50000
    HastagPackage = 500

    Transmissionrange, P_main, P_side = diskrange(systemvariables)
    print('disk: %.5f ' % (Transmissionrange))

    # initializing graph class
    ex = GraphCreation2D(x1=-10, x2=10, y1=-10, y2=10,
                         intensity=Transmissionvariables['intensity'])

    # start simulations
    nsim = 1
    for i in range(nsim):

        # create graph
        Graph, Location, N, D = ex(Tr=Transmissionrange, SD=SourceDestination,
                                   Tv=Transmissionvariables)

        # x-coordinates for all devices
        x_cord = np.array([Location[i][0] for i in range(N)])
        # y-coordinates for all devices
        y_cord = np.array([Location[i][1] for i in range(N)])
        # probability of a device interfering
        P_ = ((systemvariables['width']/(2*np.pi)) * p_trans
              * ((systemvariables['width']/(2*np.pi)) * P_main
                 + (1-(systemvariables['width']/(2*np.pi))) * P_side))

        mathModelInt = {}
        for key in D:
            # for each connected pair:
            #   one over the probability of interference being above limit.
            # Used as weight
            mathModelInt[key] = 1/sc.special.erfc(
                    np.sqrt((P_*(np.pi**3)*intensity**2/2) /
                            (_interferencelimit(
                                    systemvariables['R'], P_main, D[key]) * 2))
                    )

        # give each edge a weight, cdf
        nx.set_edge_attributes(Graph, mathModelInt, 'cdf')

        print('Connected Graph?: ' + repr(nx.is_connected(Graph)) + repr(i))

        if nx.has_path(Graph, 'S', 'D'):
            # create paths by weight, with build in algorithm
            path0 = nx.shortest_path(Graph, source='S', target='D')
            path1 = nx.shortest_path(Graph, source='S', target='D',
                                     weight='cdf')

            route_edges0 = [(path0[n], path0[n+1])
                            for n in range(len(path0)-1)]
            route_edges1 = [(path1[n], path1[n+1])
                            for n in range(len(path1)-1)]

            # conduction of transmissions, loops until all packets is
            # transmittet or max iter is achieved
            #
            # Special reuse is utilised.
            #
            # Returns number of transmitted packets and the number of
            # time slots used.
            p_received0, t_used0 = system_transmission(
                         P_main, P_side, P_,
                         p_trans, systemvariables, Location, D,
                         path0, maxiter, HastagPackage, Transmissionvariables)
            p_received1, t_used1 = system_transmission(
                         P_main, P_side, P_,
                         p_trans, systemvariables, Location, D,
                         path1, maxiter, HastagPackage, Transmissionvariables)

            # renameing (?)
            sp_received0 = p_received0
            sp_received1 = p_received1
            st_used0 = t_used0
            st_used1 = t_used1

            # Determine resulting Throughput and latency
            # the denominator transform into [packet/s] as the model from
            # [packet/timeslot]
            spath0ptr = (p_received0/t_used0)/(1/Transmissionvariables['R']+1)
            spath1ptr = (p_received1/t_used1)/(1/Transmissionvariables['R']+1)
            # Same transformation
            latency0 = (t_used0/p_received0)*(1/Transmissionvariables['R']+1)
            latency1 = (t_used1/p_received1)*(1/Transmissionvariables['R']+1)

            # summation for determination of average
            throughput_sum0 += spath0ptr
            latency_sum0 += latency0
            throughput_sum1 += spath1ptr
            latency_sum1 += latency1

            print('time slot1: %.5f, timeslot2: %.5f, thp: %.5f, %.5f'
                  % (st_used0, st_used1, spath0ptr, spath1ptr))

            # Determine corresponding expected results from model for
            # comparison
            Pt_0 = np.zeros(len(path0)-1)
            var0 = 0
            # probability of transmission for each hop
            for i in range(len(path0)-1):
                Pt_0[i] = float(
                 sc.special.erfc(
                         np.sqrt(
                                 (P_*(np.pi**3)*intensity**2/2)
                                 / (2*((P_main*(D.get(
                                    (path0[i], path0[i+1]), D.get(
                                            (path0[i+1], path0[i]))))**(-4))
                                    / (2**(Transmissionvariables['R'])-1))-1))
                                ))
                var0 += 1/Pt_0[i]           # sum of invers probabilities

            Pt_1 = np.zeros(len(path1)-1)
            var1 = 0
            for i in range(len(path1)-1):
                Pt_1[i] = float(
                 sc.special.erfc(
                         np.sqrt(
                                 (P_*(np.pi**3)*intensity**2/2)
                                 / (2*((P_main*(D.get(
                                    (path1[i], path1[i+1]), D.get(
                                            (path1[i+1], path1[i]))))**(-4))
                                    / (2**(Transmissionvariables['R'])-1))-1))
                                ))
                var1 += 1/Pt_1[i]

            # Expected throughput with spetial reuse
            ET0 = ((len(path0)-1)/2) * (1/((1/Transmissionvariables['R'] + 1)
                                        * var0))
            ET1 = ((len(path1)-1)/2) * (1/((1/Transmissionvariables['R'] + 1)
                                        * var1))
            # Expected throughput without spetial
            ET3 = (1/((1/Transmissionvariables['R']+1)*var0))
            ET4 = (1/((1/Transmissionvariables['R']+1)*var1))
            # Expected  Latency
            EL0 = ((1/Transmissionvariables['R']+1)*var0)
            EL1 = ((1/Transmissionvariables['R']+1)*var1)

            # comparison with simulation, summation is done to determine
            # the average

            # path0
            # abselute difference Throughput
            abs0 += np.abs(spath0ptr-ET0)
            # percentage deviation Throughput
            pro0 += ((np.abs(spath0ptr-ET0))/ET0)*100
            # above lower limet? (without spe.)
            if spath0ptr >= ET4:
                above0 += 1

            # abselute difference Latency
            l_abs0 += np.abs(latency0-EL0)
            # percentage deviation Latency
            l_pro0 += ((np.abs(latency0-EL0))/EL0)*100

            # path1
            abs1 += np.abs(spath1ptr-ET1)
            pro1 += ((np.abs(spath1ptr-ET1))/ET1)*100
            if spath1ptr >= ET4:
                above1 += 1

            l_abs1 += np.abs(latency1-EL1)
            l_pro1 += ((np.abs(latency1-EL1))/EL1)*100

            # summation for determine average
            ET0_sum += ET0
            ET1_sum += ET1
            EL0_sum += EL0
            EL1_sum += EL1

            # connectivity count
            if nx.is_connected(Graph) is True:
                connection_sum += 1

        else:
            print("There isn't a path between Source and Destination")

    # array with average results
    save0 = np.matrix([[abs0/nsim,
                        pro0/nsim,
                        throughput_sum0/nsim,
                        ET0_sum/nsim,
                        l_abs0/nsim,
                        l_pro0/nsim,
                        latency_sum0/nsim,
                        EL0_sum/nsim,
                        above0,
                        connection_sum]])
    save1 = np.matrix([[abs1/nsim,
                        pro1/nsim,
                        throughput_sum1/nsim,
                        ET1_sum/nsim,
                        l_abs1/nsim,
                        l_pro1/nsim,
                        latency_sum1/nsim,
                        EL1_sum/nsim,
                        above1,
                        connection_sum]])

    # save results in txt file, remember to change name!!!!
#    np.savetxt('data/R2_I8_0.out', save0, delimiter=',')
#    np.savetxt('data/R2_I8_1.out', save1, delimiter=',')

# prints garph of the last simulation for illustration.
colour = np.array(['C0' for i in range(len(Graph)-2)])
colour = np.hstack((colour, ['C1', 'C1']))
nx.draw_networkx_nodes(Graph, Location, with_labels=False, node_size=20,
                       node_color=colour, edge_color='slategrey')
nx.draw_networkx_edges(Graph, Location, route_edges0, edge_color='C2',
                       width=4)
nx.draw_networkx_edges(Graph, Location, route_edges1, edge_color='C3',
                       width=4)
