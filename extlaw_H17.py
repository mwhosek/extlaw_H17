import numpy as np
from scipy import interpolate
import pylab as py


def extinction(AKs, wavelength):
    """
    Returns the total extinction A_lambda and corresponding error (1-sigma) at the 
    specified wavelength and overall A_Ks according to the Hosek+17 extinction law. 
    This law is derived using the combined Wd1 main-squence and Arches field RC star samples, 
    assuming a Wd1 cluster age of 5 Myr. Changing the cluster age between 4 Myr - 7 Myr has no  
    effect on the law.

    Parameters:
    ----------
    AKs: float
        The overall extinction at the Ks-band (central wavelength: 2.14 microns).
        Units are magnitudes

    wavelength: float or array
        Wavelength or array of wavelengths to calculate A_lambda at. Units are
        microns. Must be between 0.8059 -- 2.2 microns, otherwise law is not defined 

    Output:
    -------
    A_lambda: array with total extinction at input wavelength(s). If single wavelength given, then
    output is a float. If array of wavelegths given, then output is a numpy array of 
    A_lambda values corresponding to the input.  
    
    err: array with the 1-sigma error. This error is symmetric about A_lambda

    Example:
    -------
    To get total extinction at 1.25 microns at overall A_Ks = 0.5 mags, with errors:
    
    >import extlaw_H17
    >result, err = extlaw_H17.extinction(0.5, 1.25) 
    
    > result[0] # extinction at 1.25 microns
    > result[0] + err[0] # extinction + 1-sigma err
    > result[0] - err[0] # extinction - 1-sigma err

    """
    # If input wavelength float, turn into array
    if type(wavelength) is float:
        wavelength = np.array([wavelength])
        
    # Check wavelegth range
    if ((min(wavelength) < 0.8) | (max(wavelength) > 2.2)):
        msg = 'Extinction law not defined at wavelength. Please select value between 0.8 - 2.2 microns'
        raise ValueError(msg)

    # Define extinction law (A_lambda / A_Ks) according to Hosek+17
    wave_law = np.array([0.8059, 0.962, 1.25, 1.53, 2.14, 3.545])
    A_AKs_law = np.array([9.66, 6.29, 3.56, 2.33, 1.0, 0.50])

    # Define the 1-sigma envelope for extinction law, according to Hosek+17
    A_AKs_low = np.array([9.34, 6.10, 3.46, 2.27, 1.0, 0.53])
    A_AKs_high = np.array([9.98, 6.48, 3.66, 2.39, 1.0, 0.47])

    # Interpolate over the curve with cubic spline interpolation
    spline_interp = interpolate.splrep(wave_law, A_AKs_law, k=3, s=0)
    spline_interp_high = interpolate.splrep(wave_law, A_AKs_high, k=3, s=0)
    spline_interp_low = interpolate.splrep(wave_law, A_AKs_low, k=3, s=0)

    # Evaluate law at input wavelength(s)
    A_AKs_at_wave = interpolate.splev(wavelength, spline_interp)
    A_AKs_at_wave_high = interpolate.splev(wavelength, spline_interp_high)
    A_AKs_at_wave_low = interpolate.splev(wavelength, spline_interp_low)

    # Multiply A_lambda / A_Ks * A_Ks to produce total extinction
    A_lambda = A_AKs_at_wave * AKs
    A_lambda_high = A_AKs_at_wave_high * AKs
    A_lambda_low = A_AKs_at_wave_low * AKs
    
    # Return A_lambda as well as error. Error is the average of the
    # positive error and negative error. 
    err_high = A_lambda_high - A_lambda
    err_low = A_lambda - A_lambda_low
    err = np.mean([err_high, err_low], axis=0)

    return A_lambda, err

def plot_extinction_law():
    """
    Helper function to plot the extinction law (with errors) over the
    wavelength range (0.8059 -- 2.2 microns). 

    Saves plot in current working directory as 'extlaw_H17.png'
    """
    # Define wavelengths
    wave = np.arange(0.8, 2.2, 0.1)

    # Define extinction law and errors. Matching derivation in paper,
    # A_Ks = 1
    law, err = extinction(1.0, wave)
    max_law = law + err
    min_law = law - err

    # Plot law
    py.figure(figsize=(10,10))
    py.plot(1.0 / wave, law, 'r-', linewidth=2, label='Wd1+RC')
    py.fill_between(1.0/wave, max_law, min_law, alpha = 0.3, facecolor = 'red',
                        label='1-sigma')
    py.xlabel(r'1 / $\lambda$ ($\mu$m$^{-1}$)', fontsize=28)
    py.ylabel(r'A$_{\lambda}$ / A$_{Ks}$', fontsize=28)
    py.text(0.46, 0.85, 'Ks', fontweight='bold')
    py.text(0.61, 0.85, 'F160W', fontweight='bold')
    py.text(0.75, 0.85, 'F125W', fontweight='bold')
    py.text(1.12, 0.85, 'F814W', fontweight='bold')
    py.text(1.01, 0.85, 'Y', fontweight='bold')
    py.legend(loc=2)
    py.savefig('extlaw_H17.png')

    return
