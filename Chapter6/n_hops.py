# -*- coding: utf-8 -*-
"""
This file is created by group C4-202 Mattek6

The following script is used for analysis the expected throughput with/witout
spatial reuse and expected latency for n-hops connections and with varying
densities.

To see the results uncomment the different sections one at time.
"""

import os
import numpy as np
from scipy.special import erfc
import matplotlib.pyplot as plt

# =============================================================================
# Function for calculating the throughput
# =============================================================================
def _interferencelimit(R, power=1, distance=1):
    """
    This function calculate the interference limit theta
    
    Input
        R:        Transmission rate R
        power:    Power P of the devices
        distance: The distance d between the devices
    
    Output
        Theta:    The interference limit   
    """
    signal = power * (distance**-4)
    return signal/(2**R-1)  - 1.

def TransmissionsSystem(R, Data, Theta, ProcessingTime, efficiency,
                        Lambda, distance, NumberHops):
    """
    This function calculate the latency, transmission success probability,
    throughput for one-hop, throughput for n-hops and throughput with spatial
    reuse.
    
    Input
        R:              Transmission rate R
        Data:           The size D of a package in bits
        Theta:          The beam width theta in degree
        ProcessingTime: The processing time Lp 
        efficiency:     Giving by nu and is use for the pratical beam width
        Lambda:         The denisty lambda of the WMN
        distance:       The distance d between the devices
        NumberHops:     The amount of hops n in the WMN
    
    Output
        Systeminfo:     A dictionary with the values for latency, transmission 
                        success probability, throughput for one-hop, throughput 
                        for n-hops and throughput with spatial reuse
    """
    TimeSlots = (Data/R) + ProcessingTime
    G_max, G_min = gaincalculation(np.radians(Theta), efficiency)

    p_theta = np.radians(Theta)/(2*np.pi)
    P = p_theta * (p_theta * G_max + (1-p_theta) * G_min)

    constant = (np.pi**3 * Lambda**2 * P)/2.

    limit = np.sqrt(constant/((_interferencelimit(R, G_max,
                                                 float(distance/NumberHops)))*2))
    
    ThroughputNspacial = (1/2) * erfc(limit) * Data * (1/TimeSlots)
    ThroughputN = (1/NumberHops) * erfc(limit) * Data * (1/TimeSlots)
    Throughput = erfc(limit) * Data * (1/TimeSlots)
    Latency = 1./erfc(limit) * NumberHops * TimeSlots
    TSP = erfc(limit)

    Systeminfo = {'Latency': Latency, 'TSP': TSP, 'Throughput': Throughput,
                  'ThroughputN': ThroughputN,
                  'ThroughputNs': ThroughputNspacial}
    return Systeminfo

def gaincalculation(width, efficiency):
    """
    This function calculate the gain g1 and g2 which are used for main and side
    lobes of a device.
    
    Input
        width:      The beam width in degree
        efficiency: Giving by nu and is use for the pratical beam width
    
    Output
        g1, g2:     The gains for main and side lobes
    """
    g1 = 2./(1. - np.cos(width/(efficiency*2.)))
    Delta = 2./(1 - np.cos(width/2.))
    g2 = (Delta - g1)/(Delta - 1.)
    return g1, g2

# =============================================================================
# Calculating the latency for varying parameters
# =============================================================================
if __name__ == '__main__':
    # import time

    """ initial and variables """
    R = 2               # From transfer plot R = 2
    Data = 1            # Amount data in bits
    Theta = 30          # Beam width in degree
    ProcessingTime = 1  # Processing latency Lp
    Efficiency = 0.8    

    Lambda = 1          # Density
    Distance = 1        # Distance between devices
    Lin = 100000
    
    """ Calculate n hops for throughput with density equal 1, 2 and 3 """
    LinNumberHops = np.linspace(1, 100, Lin)  # Number of hops from 1 to 100 hops

    Throughput1_lambda1 = np.zeros(Lin)
    Throughput1_lambda2 = np.zeros(Lin)
    Throughput1_lambda3 = np.zeros(Lin)
    
    Throughput2_lambda1 = np.zeros(Lin)
    Throughput2_lambda2 = np.zeros(Lin)
    Throughput2_lambda3 = np.zeros(Lin)

    Latency1_lambda1 = np.zeros(Lin)
    Latency1_lambda2 = np.zeros(Lin)
    Latency1_lambda3 = np.zeros(Lin)

    for i in range(Lin):
        System1_lambda1 = TransmissionsSystem(R, Data, Theta, ProcessingTime,
                                              Efficiency, 1, Distance,
                                              LinNumberHops[i])
        System1_lambda2 = TransmissionsSystem(R, Data, Theta, ProcessingTime,
                                              Efficiency, 2, Distance,
                                              LinNumberHops[i])
        System1_lambda3 = TransmissionsSystem(R, Data, Theta, ProcessingTime,
                                              Efficiency, 3, Distance,
                                              LinNumberHops[i])

        Throughput1_lambda1[i] = System1_lambda1['ThroughputNs']
        Throughput1_lambda2[i] = System1_lambda2['ThroughputNs']
        Throughput1_lambda3[i] = System1_lambda3['ThroughputNs']
        
        Throughput2_lambda1[i] = System1_lambda1['ThroughputN']
        Throughput2_lambda2[i] = System1_lambda2['ThroughputN']
        Throughput2_lambda3[i] = System1_lambda3['ThroughputN']

        Latency1_lambda1[i] = System1_lambda1['Latency']
        Latency1_lambda2[i] = System1_lambda2['Latency']
        Latency1_lambda3[i] = System1_lambda3['Latency']


    """ Expected throughput with spatial reuse """
#    plt.title('Expected Throughput', size=14)
#    plt.text(70, 0.10, r'$Rate=%s, \ \phi = %s$' % (R, Theta))
#    #plt.text(14, 0.09, r'$Rate=%s, \ \phi = %s$' % (R, Theta)) # Scale version text
#    plt.plot(LinNumberHops, Throughput1_lambda1,
#             label="Throughput with $\lambda$=1")
#    plt.plot(LinNumberHops, Throughput1_lambda2,
#             label="Throughput with $\lambda$=2")
#    plt.plot(LinNumberHops, Throughput1_lambda3,
#             label="Throughput with $\lambda$=3")
#    #plt.axis([0,20,0,0.35])                                     # Scale version
#    plt.legend()
#    plt.xlabel('Hops', size=14)
#    #plt.ylabel('Throughput [bits/s]', size=14)
#    plt.savefig("figures/throughput_nhops.PNG")
#    #plt.savefig("figures/throughput_nhops_scale.PNG")           # Scale version
#    plt.show()

    """ Expected throughput without spatial reuse """
#    plt.title('Expected Throughput', size=14)
#    plt.text(70, 0.27, r'$Rate=%s, \ \phi = %s$' % (R, Theta))
#    #plt.text(14, 0.28, r'$Rate=%s, \ \phi = %s$' % (R, Theta)) # Scale version text
#    plt.plot(LinNumberHops, Throughput2_lambda1,
#             label="Throughput with $\lambda$=1")
#    plt.plot(LinNumberHops, Throughput2_lambda2,
#             label="Throughput with $\lambda$=2")
#    plt.plot(LinNumberHops, Throughput2_lambda3,
#             label="Throughput with $\lambda$=3")
#    #plt.axis([0,20,0,0.40])
#    plt.legend()
#    plt.xlabel('Hops', size=14)
#    plt.ylabel('Throughput [bits/s]', size=14)
#    plt.savefig("figures/throughput_notspatial.PNG")
#    #plt.savefig("figures/throughput_notspatial_scale.PNG")  # Scale version
#    plt.show()    
    
    """ Expected latency """   
#    plt.title('Expected Latency', size=14)
#    plt.text(2.4,25, r'$Rate=%s, \ \phi = %s$' % (R,Theta)) # Scale version text
#    #plt.text(5, 112, r'$Rate=%s, \ \phi = %s$' % (R, Theta))
#    plt.plot(LinNumberHops, Latency1_lambda1, label="Latency with $\lambda$=1")
#    plt.plot(LinNumberHops, Latency1_lambda2, label="Latency with $\lambda$=2")
#    plt.plot(LinNumberHops, Latency1_lambda3, label="Latency with $\lambda$=3")
#    plt.axis([0,20,0,35])                                     # Scale version  
#    plt.legend()
#    plt.xlabel('Hops', size=14)
#    plt.ylabel('Latency [s]', size=14)
#    #plt.savefig("figures/latency_nhops.PNG")
#    plt.savefig("figures/latency_nhops_scale.PNG")            # Scale version
#    plt.show()
