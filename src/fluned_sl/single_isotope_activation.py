"""
this module contains some general function to calculate the decay and activation
of a single isotope in a circuit node
"""
import numpy as np


def outlet_activity(a_in_vol,rr_vol,lambda_decay,res_time):
    """
    this function calculates the activity of a node at the outlet
    """

    a_out_vol = (   rr_vol*(1-np.exp(-lambda_decay*res_time))
                  + a_in_vol*np.exp(-lambda_decay*res_time)
                )

    return a_out_vol

def rr_outlet_activity(rr_vol,lambda_decay,res_time):
    """
    this function calculates the contribution to the outlet activity of a node
    of the function value at a certain instant t. It is the given by the
    common outlet formula differentiated by the time
    """

    a_out_vol = lambda_decay*rr_vol*np.exp(-lambda_decay*res_time)

    return a_out_vol

def average_activity(a_in_vol,rr_vol,lambda_decay,res_time):
    """
    this function calculates the average activity of a node
    """

    if res_time == 0:
        a_avg_vol = a_in_vol + rr_vol
        return a_avg_vol

    a_avg_vol = (   1/res_time *
             (
             rr_vol*res_time +
           + (1/lambda_decay)*(rr_vol-a_in_vol)*(np.exp(-lambda_decay*res_time))
           - (1/lambda_decay)*(rr_vol-a_in_vol)
             )
                )


    return a_avg_vol

def rr_average_activity(rr_vol,lambda_decay,res_time):
    """
    this function calculates the contribution to the average activity of a node
    of the function value at a certain instant t. It is the given by the
    common average formula differentiated by the time
    """



    a_avg_vol = ((rr_vol/(lambda_decay*(res_time**2))) *
                    (1-np.exp(-lambda_decay*res_time)*
                    (1+lambda_decay*res_time)))


    return a_avg_vol

def res_time_from_reduction_rate(red_rate,lambda_decay):
    """
    this function calculates the residence time using the lambda decay and
    the reduction rate
    """

    res_time = - np.log(red_rate)/lambda_decay

    return res_time



