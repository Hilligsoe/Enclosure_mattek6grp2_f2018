# -*- coding: utf-8 -*-
"""
Created on Thu May 03 15:56:31 2018

@author: Frederik Vardinghus
"""

import time
import numpy as np
import networkx as nx
from main import main
from path_import import diskrange, GraphCreation2D, sufficient_range
from numba import jit
import warnings
warnings.filterwarnings("ignore")

@jit
def _interferencelimit(R, power=1, distance=1):
    signal = power * (distance**-4)
    return signal/(2**R-1) - 1.


def transmissions_success(
                     x_cord, y_cord, trans, receiv, P_main, P_side, P_prob,
                     density, width, p_trans, systemvariables, Location, D):
    trans_x, trans_y = Location[trans]
    receiv_x, receiv_y = Location[receiv]
    interfer, propinterfer = main(x_cord, y_cord, trans_x, trans_y, receiv_x,
                                  receiv_y, P_main, P_side, P_prob, density,
                                  width, p_trans)
    limit = _interferencelimit(systemvariables['R'],
                               systemvariables['power']*P_main,
                               D.get((trans, receiv), D.get((receiv, trans))))
    
    if interfer < limit:
        return 1
    else:
        return 0



def system_transmission(
                     P_main, P_side, P_,
                     p_trans, systemvariables, Location, D,
                     path, maxiter, HastagPackage, Transmissionvariables):
    """
    Function that conduct transmission of all packets throgh a given path
    in a given graph, loops until all packets is at the receiver.
    Returns number of packets at receiver and number of used time slots
    
    """
    datapath = np.zeros(len(path))
    datapath[0] = HastagPackage
#    t = time.time()
    for i in range(maxiter): # time slots
        tm = np.greater(datapath[:-1], 0)   # transmitters  (i=0 first entry true)
        rev_tm = tm[::-1]                   #               (i=0 last entry true)
        for k in range(len(rev_tm)-1):      # len is number packets -1
            if rev_tm[k] == True and rev_tm[k+1] == True:
                rev_tm[k+1] = False         # two adjencent node can't both transmit
        nodestransmitting = np.where(tm == True)[0] # array with index of transmitters
        sfT = np.zeros(len(nodestransmitting)) 
        for j in range(len(nodestransmitting)): # for each transmitting node
            trans = path[nodestransmitting[j]]      # transmitter
            receiv = path[nodestransmitting[j]+1]   # receiver
            indsys = {} 
            for key in Location: # a number assignd to each device
                if key != trans or key != receiv:
                    indsys[key] = Location[key] # kopi of Location dictionary without trans and receiver

            cord = list(indsys.values()) # list with tuples with coordinates to ever device in indsys
            x_cord = np.array([cord[h][0] for h in range(len(cord))])
            y_cord = np.array([cord[h][1] for h in range(len(cord))])

            sfT[j] = transmissions_success(x_cord, y_cord, trans, receiv,
                                           P_main, P_side, P_,
                                           Transmissionvariables['intensity'],
                                           systemvariables['width'],
                                           p_trans,
                                           systemvariables, Location, D)    
                    # returns 1 for transmission and 0 for fail to each 
                    # transmitting node in the timeslot

        datapath[nodestransmitting[np.flatnonzero(sfT)]] -= 1 # count down the number packets, for the indices with success transmission
        datapath[nodestransmitting[np.flatnonzero(sfT)]+1] += 1 # send the packets to the next device
        if datapath[-1] == HastagPackage: # done
#            print('transmission sim time: ' + repr(time.time()-t))
            return datapath[-1], i+1
        
#    print('transmission sim time: ' + repr(time.time()-t))
    return datapath[-1], i+1 # number of packets at receiver, number of time slots used


def f(intensity):
    return intensity


def Rate(n, P_main):
    return np.log2(P_main*(np.sqrt((11*np.log(n))/(n*np.pi)))**(-4) + 1)


#if __name__ == '__main__':
##    np.random.seed(3)
#
#    systemvariables = {'width': np.radians(30), 'efficiency': 0.8, 'R': 4,
#                       'power': 1, 'noise': 1}
#    SourceDestination = {'S': (-5, 0), 'D': (5, 0)}
#    Transmissionvariables = {'R': systemvariables['R'], 'packagesize': 1,
#                             'Processinglag': 1, 'intensity': 5}
#    p_trans = 0.5
#    maxiter = 5000
#    HastagPackage = 500
#    # Transmissions setup and graph
#    Transmissionrange, P_main, P_side = diskrange(systemvariables)
#
#    ex = GraphCreation2D(x1=-5, x2=5, y1=-5, y2=5,
#                         intensity=Transmissionvariables['intensity'])
#    spath0ptr = 0  
#    spath1ptr = 0
#    sp_received0 = 0
#    sp_received1 = 0
#    st_used0 = 0
#    st_used1 = 0
#    succesfullt = 0
#    nsim = 100
#    for i in range(nsim):
#        Graph, Location, N, D = ex(Tr=Transmissionrange, SD=SourceDestination,
#                                   Tv=Transmissionvariables)
#
#        P_ = ((systemvariables['width']/(2*np.pi)) * p_trans
#              * ((systemvariables['width']/(2*np.pi)) * P_main
#                 + (1-(systemvariables['width']/(2*np.pi))) * P_side))
#        s_range = sufficient_range(Transmissionvariables['intensity'], f)
#
#        print('Connected Graph?: ' + repr(nx.is_connected(Graph)) + repr(i))
#        if nx.has_path(Graph, 'S', 'D'):
#            path0 = nx.shortest_path(Graph, source='S', target='D')
#            path1 = nx.shortest_path(Graph, source='S', target='D',
#                                     weight='distance')
#            path2 = nx.shortest_path(Graph, source='S', target='D',
#                                     weight='d4')
#            path3 = nx.shortest_path(Graph, source='S', target='D',
#                                     weight=('d4', 'hop_limit'))
#
#            route_edges0 = [(path0[n], path0[n+1]) for n in range(len(path0)-1)]
#            route_edges1 = [(path1[n], path1[n+1]) for n in range(len(path1)-1)]
#            route_edges2 = [(path2[n], path2[n+1]) for n in range(len(path2)-1)]
#            route_edges3 = [(path3[n], path3[n+1]) for n in range(len(path3)-1)]
#
#            # HACK: This.
#            p_received0, t_used0 = system_transmission(
#                         P_main, P_side, P_,
#                         p_trans, systemvariables, Location, D,
#                         path0, maxiter, HastagPackage, Transmissionvariables)
#            p_received1, t_used1 = system_transmission(
#                         P_main, P_side, P_,
#                         p_trans, systemvariables, Location, D,
#                         path2, maxiter, HastagPackage, Transmissionvariables)
#            sp_received0 += p_received0
#            sp_received1 += p_received1
#            st_used0 += t_used0
#            st_used1 += t_used1
#            spath0ptr += p_received0/t_used0
#            spath1ptr += p_received1/t_used1
#            succesfullt += 1
#            save = np.array([sp_received0, st_used0, spath0ptr, sp_received1,
#                             st_used1, spath1ptr, succesfullt, nsim])
#            np.save('simulation_data_r4_l5_test13', save)
#
#        else:
#            print("There isn't a path between Source and Destination")
