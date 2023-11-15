# Dryden Gust Model 
# Mid size modelling
# Date: 10th Nov 2023
#----------------------------------------------------------#
# Module import
import scipy.io as sio
import numpy as np
import math
from scipy.signal import butter, lfilter, freqz, firwin
from scipy import signal
from matplotlib import pyplot as plt
import csv
from matplotlib.colors import Normalize


# Low altitude model
# Transfer function for along-wind
def u_transfer_function(height, airspeed):
    # turbulence level defines value of wind speed in knots at 20 feet
    # turbulence_level = 15 * 0.514444 
    # convert speed from knots to meters per second
    turbulence_level = 15 
    Lu = height / ((0.177 + 0.000823*height) ** (1.2))
    Lu_V = Lu / airspeed
    # length_u = 1750
    sigma_w = 0.1 * turbulence_level 
    sigma_u = sigma_w / ((0.177 + 0.000823*height) ** (0.4))
    num_u = [(2 * (sigma_u ** 2) * Lu_V) / (math.pi)]
    den_u = [(Lu_V ** 2), 1]
    phi_u = signal.TransferFunction(num_u, den_u)
    return phi_u


# Transfer function for vertical-wind
def v_transfer_function(height, airspeed):
    # turbulence level defines value of wind speed in knots at 20 feet
    # turbulence_level = 15 * 0.514444 
    # convert speed from knots to meters per second
    turbulence_level = 15 
    length_v = height / ((0.177 + 0.000823*height)**(1.2))
    # length_v = 1750
    sigma_w = 0.1 * turbulence_level 
    sigma_v = sigma_w / ((0.177 + 0.000823*height) ** (0.4))
    b = sigma_v * (math.sqrt((length_v) / (math.pi * airspeed)))
    Lv_V = length_v / airspeed
    num_v = [(math.sqrt(3) * Lv_V * b), b]
    den_v = [(Lv_V**2), 2*Lv_V, 1]
    phi_v = signal.TransferFunction(num_v, den_v)
    return phi_v

def w_transfer_fucntion(height, airspeed):
    # turbulence level defines value of wind speed in knots at 20feet
    # turbulence_level = 15 * 0.514444 
    # convert speed from knots to meters per second
    turbulence_level = 15
    length_w = height
    sigma_w = 0.1 * turbulence_level
    c = sigma_w * (math.sqrt((length_w) / (math.pi * airspeed)))
    Lw_V = length_w / airspeed
    num_w = [(math.sqrt(3) * Lw_V * c), c]
    den_w = [(Lw_V ** 2), 2 * Lw_V, 1]
    H_v = signal.TransferFunction(num_w, den_w)
    return H_v

# dryden wind model forming
def dryden_wind_velocities(time, height, airspeed):
    mean = 0
    std = 1

    t_p = np.linspace(0, time, 1000)
    num_sample = 1000

    # random number seed use the same as from the simulink in matlab
    # ug
    np.random.seed(23341)
    sample_1 = 10 * np.random.normal(mean, std, size = num_sample)
    # vg
    np.random.seed(23342)
    sample_2 = 10 * np.random.normal(mean, std, size = num_sample)
    # wg
    np.random.seed(23343)
    sample_3 = 10 * np.random.normal(mean, std, size = num_sample)

    # transfer velocities
    tf_u = u_transfer_function(height, airspeed)
    tf_v = v_transfer_function(height, airspeed)
    tf_w = w_transfer_fucntion(height, airspeed)


    # compute the response to transfer function
    tout1, yo1, xo1 = signal.lsim(tf_u, sample_1, t_p)
    yo1_f = [i * 0.305 for i in yo1]
    tout2, yo2, xo2 = signal.lsim(tf_v, sample_2, t_p)
    yo2_f = [i * 0.305 for i in yo2]
    tout3, yo3, xo3 = signal.lsim(tf_w, sample_3, t_p)
    yo3_f = [i * 0.305 for i in yo3]

    # plot the along-direction wind velocities generation
    plt.figure(1)
    plt.plot(t_p, yo1_f, 'b')
    plt.ylabel('along-wind in m/s (P)')
    plt.xlabel('time in s')
    plt.grid(True)
    plt.savefig('./along-wind.jpg')

    # plot the crossing-direction wind velocities generation
    plt.figure(2)
    plt.plot(t_p, yo2_f, 'r')
    plt.ylabel('crossing-wind in m/s (P)')
    plt.xlabel('time in s')
    plt.grid(True)
    plt.savefig('./crossing-wind.jpg')

    # plot the vertical-direction wind velocities generation
    plt.figure(3)
    plt.plot(t_p, yo3_f, 'g')
    plt.ylabel('along-wind in m/s (P)')
    plt.xlabel('time in s')
    plt.grid(True)
    plt.savefig('./vertical-wind.jpg')

    # show all plots
    plt.show()



# main program
def main():
    t = int(input("Enter the time of wind generation (Unit: s): "))
    h = float(input("Enter the altitude of flight (Unit: m): "))
    V = float(input("Enter the cruising speed (Unit: m/s): "))

    # using the dryden model
    dryden_wind_velocities(t, h, V)

main()










