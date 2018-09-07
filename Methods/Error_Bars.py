# This script gets the error bars for the individual points when measuring the differential scintillation
import numpy as np
import math

def ipartial(real, imaginary):
    return (1/(real*(1 + (imaginary/real)**2)))**2

def rpartial(real, imaginary):
    return (1/(real**2*(1 +(imaginary/real)**2)))**2

def GetPlot(cross, xlim):
    """Return the plot  errorbars for the associated points as a one-dimensional numpy array

    Cross: The cross secondary spectrum, already binned. First dimension assumed to be 
    doppler frequency, second dimension is time delay. Only positive time delays should be
    plotted.
    """
    
    copy = RemoveNans(cross)
    mean = np.nanmean(copy, axis=0)
    real = np.real(mean)
    imaginary = np.imag(mean)
    realstd = np.nanstd(np.real(copy), axis=0)/np.sqrt(copy.shape[0])
    imagstd = np.nanstd(np.imag(copy), axis=0)/np.sqrt(copy.shape[0])
    std = np.sqrt(ipartial(real,imaginary)*(np.abs(imagstd))**2 + \
                      rpartial(real,imaginary)*(np.abs(realstd))**2)
    domain = np.linspace(start=-xlim, stop=xlim, num=std.shape[0])
    return domain, np.angle(mean), std

def GetPlot_Modified(cross, xlim, off_std):
    """Return a modified plot errorbars for the associated points as a one-dimensional numpy
    array. The difference in power between the parabola and outside the parabola is removed
    in the calculation of the error by dividing out the standard deviation of the off gate
    divided by route n ()
    Cross: The cross secondary spectrum, already binned. First dimension assumed to be doppler
    frequency, second dimension is time delay. Only positive time delays should be plotted"""

    copy = RemoveNans(cross) 
    mean = np.nanmean(copy, axis=0) # Average along the doppler time axis, collapsing frequency
    real = np.real(mean)
    imaginary = np.imag(mean)
    power_diff = off_std/np.sqrt(copy.shape[0]) # 
    realstd = np.nanstd(np.real(copy)/power_diff, axis=0)/np.sqrt(copy.shape[0])
    imagstd = np.nanstd(np.imag(copy)/power_diff, axis=0)/np.sqrt(copy.shape[0])
    std = np.sqrt(ipartial(real,imaginary)*(np.abs(imagstd))**2 + \
                      rpartial(real,imaginary)*(np.abs(realstd))**2)
    domain = np.linspace(start=-xlim, stop=xlim, num=std.shape[0])
    return domain, np.angle(mean), std


def RemoveNans(cross, xlim=None):
    """Remove vertical axis whose values are all nan and return new matrix"""
    # Temporary locations of all nan values in the vertical and horizatonal axis
    hnans, vnans = [], []
    for i in range(cross.shape[0]):
        if (np.isnan(cross[i,:])).all():
            hnans.append(i)
    for j in range(cross.shape[1]):
        if (np.isnan(cross[:,j])).all():
            vnans.append(j)
    # Delete nan axis
    copy = np.delete(cross, vnans, axis=1)
    copy = np.delete(copy, hnans, axis=0)
    length = math.ceil(cross.shape[1]/2)
    # Make the domain axis if the x axis limit is given
    if xlim is not None:
        domain = np.arange(-length, length) + 0.5
        domain *= (xlim/domain[-1])
        domain = np.delete(domain, vnans)
        return copy, domain
    else:
        return copy
    
    
    
