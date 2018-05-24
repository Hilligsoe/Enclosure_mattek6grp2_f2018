# -*- coding: utf-8 -*-
"""
Created on Thu May 03 15:56:31 2018

@author: MATTEK6 grp 2
"""

from path_import import _gaincalculation
from wmn_main import _interferencelimit
from main import main
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")


def transmissions_success(
                     x_cord, y_cord, distance, P_main, P_side, P_prob,
                     density, width, p_trans, systemvariables):
    trans_x, trans_y = 1, 0
    receiv_x, receiv_y = 0, 0
    interfer, propinterfer = main(x_cord, y_cord, trans_x, trans_y, receiv_x,
                                  receiv_y, P_main, P_side, P_prob, density,
                                  width, p_trans)
    limit = _interferencelimit(systemvariables['R'],
                               systemvariables['power']*P_main, distance)

    if interfer < limit:
        return 1
    else:
        return 0


def ppp(density, size):
    number = np.random.poisson(density*size**2)
    x = np.random.uniform(-size/2., size/2., number)
    y = np.random.uniform(-size/2., size/2., number)

    return x, y


if __name__ == '__main__':
    # initializing varibles for summation of results to determine averge
    abs0 = 0
    pro0 = 0
    l_abs0 = 0
    l_pro0 = 0
    throughput_sum0 = 0
    latency_sum0 = 0
    ET0_sum = 0
    EL0_sum = 0

    # set random seed (3 is used for documented results)
#    np.random.seed(3)

    # intital system settings
    size = 20
    intensity = 5
    distance = 1
    width = np.radians(30)
    efficiency = 0.8
    P_main, P_side = _gaincalculation(width, efficiency)
    systemvariables = {'width': width, 'efficiency': efficiency,
                       'R': 2, 'power': 1, 'noise': 1}
    Transmissionvariables = {'R': systemvariables['R'], 'packagesize': 1,
                             'Processinglag': 1, 'intensity': intensity}
    p_trans = 0.5
    maxiter = 50000
    HastagPackage = 100
    linLength = 200
    # initializing graph class
    intensity = np.linspace(0, 100, linLength)
    throughput = np.zeros(linLength)
    latency = np.zeros(linLength)
    # start simulations
    for d in range(linLength):
        nsim = 400
        throughput_sum1 = 0
        latency_sum1 = 0
        print(d)
        for i in range(nsim):

            # x, y cord for all devices
            x, y = ppp(intensity[d], size)
            # probability of a device interfering
            P_ = ((systemvariables['width']/(2*np.pi)) * p_trans
                  * ((systemvariables['width']/(2*np.pi)) * P_main
                     + (1-(systemvariables['width']/(2*np.pi))) * P_side))

            p_received = 0
            for j in range(maxiter):
                p_received += transmissions_success(
                         x, y, distance, P_main, P_side, P_,
                         intensity[d], width, p_trans, systemvariables)
                if p_received == HastagPackage:
                    break
            # Determine resulting Throughput and latency
            # the denominator transform into [packet/s] as the model from
            # [packet/timeslot]
            spath0ptr = (p_received/(j+1))/(1/systemvariables['R']+1)
            # Same transformation
            latency0 = ((j+1)/p_received)*(1/systemvariables['R']+1)

            # summation for determination of average
            throughput_sum0 += spath0ptr
            latency_sum0 += latency0
            throughput_sum1 += spath0ptr
            latency_sum1 += latency0

            # Determine corresponding expected results from model for
            # comparison
            # probability of transmission for each hop
#            Pt_0 = float(sc.special.erfc(
#                     np.sqrt(
#                             (P_*(np.pi**3)*intensity[d]**2/2)
#                             / (2*((P_main*(distance)**(-4))
#                                / (2**(systemvariables['R'])-1))-1))
#                            ))
#            var0 = 1/Pt_0           # sum of invers probabilities

            # Expected throughput with spetial reuse
#            ET0 = (1/((1/systemvariables['R'] + 1) * var0))
#            # Expected  Latency
#            EL0 = ((1/systemvariables['R']+1)*var0)

            # comparison with simulation, summation is done to determine
            # the average

            # path0
            # abselute difference Throughput
#            abs0 += np.abs(spath0ptr-ET0)
#            # kvocentage deviation Throughput
#            pro0 += (ET0/spath0ptr)
#
#            # abselute difference Latency
#            l_abs0 += np.abs(latency0-EL0)
#            # kvocantage deviation Latency
#            l_pro0 += EL0/latency0
#
#            # summation for determine average
#            ET0_sum += ET0
#            EL0_sum += EL0
        throughput[d] = throughput_sum1/nsim
        latency[d] = latency_sum1/nsim
        # array with average results
        save = np.matrix([[abs0/nsim,
                           pro0/nsim,
                           throughput_sum0/nsim,
                           ET0_sum/nsim,
                           l_abs0/nsim,
                           l_pro0/nsim,
                           latency_sum0/nsim,
                           EL0_sum/nsim]])

        # save results in txt file, remember to change name!!!!
    #    np.savetxt('data/1hop_R_I_x.out', save, delimiter=',')

#    plt.figure(2)
#    plt.plot(intensity, latency, label='1 hop connection')
#    plt.text(0, 1.648, r'$Rate=%s, \ \theta = %s$' % (2, 30))
#    plt.legend()
#    plt.xlabel(r'Intensity $\lambda$')
#    plt.ylabel('Latency [s]')
#    plt.title('Calculated Throughput')
#    plt.show()
#
#    plt.figure(1)
#    plt.plot(intensity, throughput, label='1 hop connection')
#    plt.text(75, 0.663, r'$Rate=%s, \ \theta = %s$' % (2, 30))
#    plt.legend()
#    plt.xlabel(r'Intensity $\lambda$')
#    plt.ylabel('Throughput [bits/s]')
#    plt.title('Calculated Latency')
#    plt.show()
