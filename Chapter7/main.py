import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


def _points_in_beam(x_0, y_0, x, y, direction, width):
    x = x - x_0
    y = y - y_0

    angle = np.angle(x+y*1j)
    negative = (angle != abs(angle))
    angle[negative] = angle[negative] + 2*np.pi

    if direction - width/2. < 0:
        x_in1 = np.where(angle < direction + width/2.)
        x_in2 = np.where(angle > 2*np.pi + direction - width/2.)
        x_in = np.sort(np.hstack((x_in1, x_in2)))
        y_in = y[x_in] + y_0
        x_in = x[x_in] + x_0

        return x_in, y_in

    elif direction + width/2. > 2*np.pi:
        x_in1 = np.where(angle < -2*np.pi + direction + width/2.)
        x_in2 = np.where(angle > direction - width/2.)
        x_in = np.sort(np.hstack((x_in1, x_in2)))
        y_in = y[x_in] + y_0
        x_in = x[x_in] + x_0

        return x_in, y_in

    else:
        x_in = np.where(
                np.logical_and(
                        np.less(
                                angle, direction+width/2.), np.greater(
                                        angle, direction-width/2.)))
        length = np.shape(x_in)[1]
        y_in = np.reshape(y[x_in] + y_0, [1, length])
        x_in = np.reshape(x[x_in] + x_0, [1, length])

        return x_in, y_in


def _direction_of_interferers(x_in, y_in, p_trans, x_0, y_0, width):
    on = np.where(np.random.binomial(1, p_trans, len(x_in[0])) == 1)
    x_on = x_in[0][on]
    y_on = y_in[0][on]

    orientation = np.random.uniform(0, 2*np.pi, len(x_on))
    angle = np.angle(x_0 - x_on + (y_0 - y_on) * 1j)
    negative = (angle != abs(angle))
    angle[negative] = angle[negative] + 2*np.pi

    orientation_check = np.logical_and(
            np.less(
                    angle, orientation + width/2.), np.greater(
                            angle, orientation - width/2.))
    x_main, y_main = x_on[orientation_check], y_on[orientation_check]
    x_side, y_side = x_on[orientation_check == 0], y_on[orientation_check == 0]

    return x_main, y_main, x_side, y_side


def _calculate_interference(
        x_receiver, y_receiver, x_main, y_main, x_side, y_side,
        P_main, density, P_side, P_prob):
    dist_main = np.sqrt((x_receiver - x_main)**2. + (y_receiver - y_main)**2.)
    dist_side = np.sqrt((x_receiver - x_side)**2. + (y_receiver - y_side)**2.)

    interference = P_main * np.sum(dist_main**(-4.)) + P_side * np.sum(
                                                            dist_side**(-4.))
    cdfI = sp.special.erfc(
            np.sqrt((P_prob*density**2.*np.pi**3/2.))/(2.*interference))

    return interference, cdfI


def _receiver_dir_dist(
        x_transmitter, y_transmitter, x_receiver, y_receiver):
    distance = np.sqrt((
            x_transmitter-x_receiver)**2 + (y_transmitter-y_receiver)**2)
    direction = np.angle(
            x_transmitter - x_receiver + (y_transmitter - y_receiver)*1j)

    if direction < 0:
        direction = direction + 2*np.pi

    return distance, direction


def main(x, y, x_transmitter, y_transmitter, x_receiver, y_receiver,
         P_main, P_side, P_prob, density, width, p_trans):

    width = width/180.*np.pi

    # =========================================================================
    # Beregn distance og retning mellem afsender og modtager
    # =========================================================================

    distance, direction = _receiver_dir_dist(
            x_transmitter, y_transmitter, x_receiver, y_receiver)

    # =========================================================================
    # Punkter i beam
    # =========================================================================

    x_in, y_in = _points_in_beam(x_receiver, y_receiver, x, y, direction, width)

    # =========================================================================
    # Sluk antenner, giv resten retning, tjek om de peger mod afsenderen
    # =========================================================================

    x_main, y_main, x_side, y_side = _direction_of_interferers(
            x_in, y_in, p_trans, x_receiver, y_receiver, width)

    # =========================================================================
    # Beregn interferens
    # =========================================================================

    interference, probability_of_interference = _calculate_interference(
            x_receiver, y_receiver, x_main, y_main, x_side, y_side, P_main,
            density, P_side, P_prob)

    # =========================================================================
    # Plots
    # =========================================================================

    # Alle punkter

#    plt.figure(1)
#    plt.scatter(x, y, c="black", marker=".")
#    plt.scatter(x_receiver, y_receiver, c="red")
#    plt.scatter(x_transmitter, y_transmitter, c="blue")
#    plt.xlim(-size/2., size/2.)
#    plt.ylim(-size/2., size/2.)

    # Alle punkter i keglen

#    plt.figure(2)
#    plt.scatter(x_in, y_in, c="black", marker=".")
#    plt.scatter(x_receiver, y_receiver, c="red")
#    plt.scatter(x_transmitter, y_transmitter, c="blue")
#    plt.xlim(-size/2., size/2.)
#    plt.ylim(-size/2., size/2.)

    # Punkter som er vendt mod afsenderen

#    plt.figure(3)
#    plt.scatter(x_main, y_main)
#    plt.scatter(x_receiver, y_receiver, c="red")
#    plt.scatter(x_transmitter, y_transmitter, c="green")
#    plt.xlim(-size/2., size/2.)
#    plt.ylim(-size/2., size/2.)

    # Punkter som ikke er vendt mod afsenderen

#    plt.scatter(x_side, y_side)
#    plt.scatter(x_receiver, y_receiver, c="red")
#    plt.scatter(x_transmitter, y_transmitter, c="green")
#    plt.xlim(-size/2., size/2.)
#    plt.ylim(-size/2., size/2.)

    return interference, probability_of_interference
