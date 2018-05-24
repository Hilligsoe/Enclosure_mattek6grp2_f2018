# -*- coding: utf-8 -*-
"""
This file is created by group C4-202 Mattek6

The following script is used for analysis the latency for varying 
transmission rate, beam width and densities.

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
    LinLength = 1000
    
    """ Calculate different transmission rates for latency """
#    LinTransferRate = np.linspace(0.01, 10, LinLength)  # Transfer rate R from 0.01 to 10
#
#    Latency1 = np.zeros(LinLength)
#    Latency2 = np.zeros(LinLength)
#
#    for i in range(LinLength):
#        System1 = TransmissionsSystem(LinTransferRate[i], Data, Theta,
#                                      ProcessingTime, Efficiency, Lambda,
#                                      Distance, NumberHops=1)
#
#        System2 = TransmissionsSystem(LinTransferRate[i], Data, Theta,
#                                      ProcessingTime, Efficiency, Lambda,
#                                      Distance, NumberHops=2)
#
#        Latency1[i] = System1['Latency']
#        Latency2[i] = System2['Latency']
#   
#    plt.title('Expected Latency', size=14)
#    plt.text(1, 24, r'$\lambda = %s, \ \phi = %s$' % (Lambda, Theta))
#    plt.plot(LinTransferRate, Latency1, label="Latency for 1 hop")
#    plt.plot(LinTransferRate, Latency2, label="Latency for 2 hops")
#    plt.legend()
#    plt.xlabel('Transmission rate [bits/s]', size=14)
#    plt.ylabel('Latency [s]', size=14)
#    plt.axis([0, 10, 0, 30])
#    plt.savefig(f"figures/latency_rate.PNG")
#    plt.show()
#    
#    """ Finding the optimal value for the transmission rate """
#    optimal = np.nanargmin(Latency1)
#    
#    plt.title('Optimal value of Expected Latency', size=14)
#    plt.text(7, 24, r'$\lambda = %s, \ \phi = %s$' % (Lambda, Theta))
#    plt.plot(LinTransferRate, Latency1, label="Latency for 1 hop")
#    plt.plot(LinTransferRate[optimal], Latency1[optimal], 'ro', label="Optimal value")
#    plt.legend()
#    plt.xlabel('Transmission rate [bits/s]', size=14)
#    plt.ylabel('Latency [s]', size=14)
#    plt.axis([0, 10, 0, 30])
#    plt.savefig(f"figures/latency_optimal.PNG")
#    plt.show()
#    
    """ Calculate different angle for Latency """
#    LinAngle = np.linspace(0, 300, LinLength)  # Beam width from 0 to 300 degree
#    
#    Latency1_angle = np.zeros(LinLength)
#    Latency2_angle = np.zeros(LinLength)
#    
#    for i in range(LinLength):
#        System1_angle = TransmissionsSystem(R, Data, LinAngle[i],
#                                      ProcessingTime, Efficiency, Lambda,
#                                      Distance, NumberHops=1)
#
#        System2_angle = TransmissionsSystem(R, Data, LinAngle[i],
#                                      ProcessingTime, Efficiency, Lambda,
#                                      Distance, NumberHops=2)
#        
#        Latency1_angle[i] = System1_angle['Latency']
#        Latency2_angle[i] = System2_angle['Latency']
#
#    plt.title('Expected Latency', size=14)
#    plt.text(220, 41, r'$Rate=%s, \ \lambda = %s$' % (R, Lambda))
#    plt.plot(LinAngle, Latency1_angle, label="Latency for 1 hop")
#    plt.plot(LinAngle, Latency2_angle, label="Latency for 2 hops")
#    plt.legend()
#    plt.xlabel(r'Angle $\phi$ [deg]', size=14)
#    plt.ylabel('Latency [s]', size=14)
#    plt.axis([0, 300, 0, 50])
#    plt.savefig("figures/latency_angle_scale.PNG")
#    plt.show()    

    """ Calculate different intensities for Latency """
#    LinDensity = np.linspace(1, 20, LinLength)  # Density from 1 to 20
#
#    Latency1_density = np.zeros(LinLength)
#    Latency2_density = np.zeros(LinLength)
#
#    for i in range(LinLength):
#        System1_density = TransmissionsSystem(R, Data, Theta, ProcessingTime,
#                                              Efficiency, LinDensity[i],
#                                              Distance, NumberHops=1)
#
#        System2_density = TransmissionsSystem(R, Data, Theta, ProcessingTime,
#                                              Efficiency, LinDensity[i],
#                                              Distance, NumberHops=2)
#
#        Latency1_density[i] = System1_density['Latency']
#        Latency2_density[i] = System2_density['Latency']   
#    
#    plt.title('Expected Latency', size=14)
#    plt.text(14, 18, r'$Rate=%s, \ \phi = %s$' % (R, Theta))
#    plt.plot(LinDensity, Latency1_density, label="Latency for 1 hop")
#    plt.plot(LinDensity, Latency2_density, label="Latency for 2 hops")
#    plt.legend()
#    plt.axis([0,20,0,100])
#    plt.xlabel(r'Intensity $\lambda$', size=14)
#    plt.ylabel('Latency [s]', size=14)
#    plt.savefig("figures/latency_density.PNG")
#    plt.show()
