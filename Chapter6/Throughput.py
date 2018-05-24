# -*- coding: utf-8 -*-
"""
This file is created by group C4-202 Mattek6

The following script is used for analysis the throughput for varying 
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
# Calculating the throughput for varying parameters
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
        
    """ Calculate different transmission rates for throughput """
#    LinTransferRate = np.linspace(0.01, 10, LinLength) # Transfer rate R from 0.01 to 10
#    
#    Throughput1 = np.zeros(LinLength)
#    Throughput2 = np.zeros(LinLength)
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
#        Throughput1[i] = System1['Throughput']
#        Throughput2[i] = System2['ThroughputNs']
#
#    plt.title('Expected Throughput', size=14)
#    plt.text(7,0.4, r'$\phi = %s, \ \lambda = %s$' % (Theta, Lambda)) # Theta = 30 
#    plt.plot(LinTransferRate, Throughput1, label="Throughput for 1 hop")
#    plt.plot(LinTransferRate, Throughput2, label="Throughput for 2 hops")
#    plt.legend()
#    plt.axis([0, 10, 0, 0.5])
#    plt.xlabel('Transmission rate [bits/s]', size=14)
#    plt.ylabel('Throughput [bits/s]', size=14)
#    plt.savefig(f"figures/throughput_rate.PNG")
#    plt.show()
#        
#    """ Finding the optimal value for the transmission rate """
#    optimal = np.nanargmax(Throughput1)
#    
#    plt.title('Optimal value of Expected Throughput', size=14)
#    plt.text(7,0.4, r'$\phi = %s, \ \lambda = %s$' % (Theta, Lambda))
#    plt.plot(LinTransferRate, Throughput1, label="Throughput for 1 hop")
#    plt.plot(LinTransferRate[optimal], Throughput1[optimal], 'ro', label="Optimal value")
#    plt.legend()
#    plt.axis([0, 10, 0, 0.5])
#    plt.xlabel('Transmission rate [bits/s]', size=14)
#    plt.ylabel('Throughput [bits/s]', size=14)
#    plt.savefig(f"figures/throughput_optimal.PNG")
#    plt.show()
#    
    """ Calculate different angle for throughput """
#    LinAngle = np.linspace(0, 300, LinLength)  # Beam width theta from 0 to 300 degree
#    
#    Throughput1_angle = np.zeros(LinLength)
#    Throughput2_angle = np.zeros(LinLength)
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
#
#
#        Throughput1_angle[i] = System1_angle['Throughput']
#        Throughput2_angle[i] = System2_angle['ThroughputNs']
#
#    plt.title('Expected Throughput', size=14)
#    plt.text(220, 0.55, r'$Rate=%s, \ \lambda = %s$' % (R, Lambda))
#    plt.plot(LinAngle, Throughput1_angle, label="Throughput for 1 hop")
#    plt.plot(LinAngle, Throughput2_angle, label="Throughput for 2 hops")
#    plt.legend()
#    plt.xlabel(r'Angle $\phi$ [deg]', size=14)
#    plt.ylabel('Throughput [bits/s]', size=14)
#    plt.savefig("figures/throughput_angle.PNG")
#    plt.show()
#    
    """ Calculate different intensities for throughput """
#    LinDensity = np.linspace(1, 20, LinLength)  # Density from 1 to 20
#
#    Throughput1_density = np.zeros(LinLength)
#    Throughput2_density = np.zeros(LinLength)
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
#        Throughput1_density[i] = System1_density['Throughput']
#        Throughput2_density[i] = System2_density['ThroughputNs']   
#    
#    plt.title('Expected Throughput', size=14)
#    plt.text(14, 0.30, r'$Rate=%s, \ \phi = %s$' % (R, Theta))
#    plt.plot(LinDensity, Throughput1_density, label="Throughput for 1 hop")
#    plt.plot(LinDensity, Throughput2_density, label="Throughput for 2 hops")
#    plt.legend()
#    plt.xlabel(r'Intensity $\lambda$', size=14)
#    plt.ylabel('Throughput [bits/s]', size=14)
#    plt.savefig("figures/throughput_density.PNG")
#    plt.show()
