###########################################################################################################
import scipy.io as sio
import numpy as np
import math
from scipy.signal import butter, lfilter, freqz, firwin
from scipy import signal
from matplotlib import pyplot as plt
import csv

# the dryden transfer function defined in MIL-HDBK-1797 and MIL-HDBK-1797B for the longtitude, lateral and vertical-wind directions:

# low altitude model (altitude < 1000ft)
# tranfer function for longtitude-wind
def u_transfer_function(height, airspeed):
    # turbulence level defines value of wind speed in knots @ 20ft: W20
    # turbulence level @ W20 is 15 knots
    W20 = 20 # cover the unit from BS to ISO: turbulence_level * 0.514444
    Lu = (height) / ((0.177 + 0.000823 * height) ** (1.2))
    sigma_w = 0.1 * W20
    Lu_V = Lu / airspeed
    sigma_u = (sigma_w) / ((0.177 + 0.000823 * height) ** (0.4))
    num_u = [(2 * (sigma_u ** 2) * Lu_V) / (math.pi)]
    den_u = [(Lu_V ** 2), 1]
    phi_u = signal.TransferFunction(num_u, den_u)
    return phi_u


# transfer function for lateral-wind
def v_transfer_function(height, airspeed):
    W20 = 20
    Lv = (height) / (2 * ((0.177 + 0.000823 * height) ** (1.2)))
    sigma_w = 0.1 * W20
    Lv_V = Lv / airspeed
    sigma_v = (sigma_w) / ((0.177 + 0.000823 * height) ** (0.4))
    a = (2 * sigma_v * Lv_V) / (math.pi)
    num_v = [(12 * a * (Lv_V ** 2)), a]
    den_v = [(16 * (Lv_V ** 4)), (8 * (Lv_V ** 2)), 1]
    phi_v = signal.TransferFunction(num_v, den_v)
    return phi_v


# transfer function for vertical-wind
def w_transfer_fucntion(height, airspeed):
    W20 = 20
    Lw = (height) / 2
    sigma_w = 0.1 * W20
    Lw_V = Lw / airspeed
    b = (2 * sigma_w * Lw_V) / (math.pi)
    num_w = [(12 * b * (Lw_V ** 2)), 0, b]
    den_w = [(16 * (Lw_V ** 4)), 0, (8 * (Lw_V ** 2)), 0, 1]
    phi_w = signal.TransferFunction(num_w, den_w)
    return phi_w


# dryden wind model by applying the transfer function
def dryden_wind_velocities_noise(time, height, airspeed):
    # height and speed from ISO unit to BS unit
    height = float(height) * 3.28084
    airspeed = float(airspeed) * 3.28084

    # generate the white noise
    mean = 0
    std = 1
    # generate the sample number list of white noise
    # create the sequence of 10000 equally spaced numerical values among 0 - 5
    numlist = np.linspace(0, time, 1000)
    num_sample = 1000

    # the random number seed used the same as from the SIMULINK in MATLAB
    # noise seed - ug
    np.random.seed(23341)
    ug_sample = 10 * np.random.normal(mean, std, size = num_sample)
    # noise seed - vg
    np.random.seed(23342)
    vg_sample = 10 * np.random.normal(mean, std, size = num_sample)
    # noise seed - wg
    np.random.seed(23343)
    wg_sample = 10 * np.random.normal(mean, std, size = num_sample)

    # generate the dryden model wind speed from three different directions
    tf_u = u_transfer_function(height, airspeed)
    tf_v = v_transfer_function(height, airspeed)
    tf_w = w_transfer_fucntion(height, airspeed)

    # compute response to transfer function
    # transfer the results unit into ISO
    tout1, y1, x1 = signal.lsim(tf_u, ug_sample, numlist)
    y1_f = [i * 0.305 for i in y1]
    tout2, y2, x2 = signal.lsim(tf_v, vg_sample, numlist)
    y2_f = [i * 0.305 for i in y2]
    tout3, y3, x3 = signal.lsim(tf_v, wg_sample, numlist)
    y3_f = [i * 0.305 for i in y3]

    # plot diagrams
    t_p = numlist
    plt.figure(1, figsize = (20, 6))
    plt.plot(t_p, y1_f, 'b', t_p, y2_f, 'r', t_p, y3_f, 'g')
    plt.legend(["longtitude-wind", "lateral-wind", "vertical-wind"])
    plt.ylabel('wind-speed in m/s (P)')
    plt.xlabel('time in s')
    plt.grid(True)
    plt.savefig('./dryden_wind_model_with_noise.jpg')

    # plot the along-direction wind velocities generation
    plt.figure(2, figsize = (20, 6))
    plt.plot(t_p, y1_f, '-b')
    plt.ylabel('along-wind in m/s (P)')
    plt.xlabel('time in s')
    plt.grid(True)
    plt.savefig('./along-wind.jpg')

    # plot the crossing-direction wind velocities generation
    plt.figure(3, figsize = (20, 6))
    plt.plot(t_p, y2_f, 'r')
    plt.ylabel('crossing-wind in m/s (P)')
    plt.xlabel('time in s')
    plt.grid(True)
    plt.savefig('./crossing-wind.jpg')

    # plot the vertical-direction wind velocities generation
    plt.figure(4, figsize = (20, 6))
    plt.plot(t_p, y3_f, 'g')
    plt.ylabel('along-wind in m/s (P)')
    plt.xlabel('time in s')
    plt.grid(True)
    plt.savefig('./vertical-wind.jpg')

    plt.show()


# main program
def main():
    t = int(input("Enter the time of wind generation (Unit: s): "))
    h = float(input("Enter the altitude of flight (Unit: m): "))
    V = float(input("Enter the cruising speed (Unit: m/s): "))

    # using the dryden model
    dryden_wind_velocities_noise(t, h, V)

main()





