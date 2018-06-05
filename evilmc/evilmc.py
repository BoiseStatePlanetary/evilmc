# -*- coding: utf-8 -*-

import numpy as np

from PyAstronomy.modelSuite.XTran.forTrans import MandelAgolLC
import transit_utils

__all__= ['evilmc']

class EVILMC(object):
    """Returns ellipsoidal variation of a slowly-rotating star induced by a 
    low-mass companion as described in Jackson et al. (2012) ApJ 750, 1 - 
    http://adsabs.harvard.edu/abs/2012ApJ...751..112J.
    """

    def __init__(self, time, params, supersample_factor=1, exp_time=0):
        """__init__ method for EVILMC.

        Args:
            time (numpy array): observational time (same units at period)

            params (:obj:`dict`): dict of floats/numpy arrays, including
                params["per"] - orbital period (any units)
                params["a"] - semi-major axis (units of stellar radius)
                params["T0"] - mid-transit time (same units as period)
                params["p"] - planet's radius (stellar radius)
                if quadratic limb-darkening,
                   params["linLimb"] - linear limb-darkening coefficient
                   params["quadLimb"] - quadratic limb-darkening coefficient
                if non-linear limb-darkening,
                    params["a*"] - limb-darkening coefficients
                params["b"] - impact parameter (stellar radius);
                    if not given, inclination angle must be
                params["q"] - planet-star mass ratio
                params["Kz"] - stellar reflex velocity (in lightspeed)
                params["Ts"] - stellar temperature (in K)
                params["Ws"] - stellar rotation vector in x, y, z
                params["beta"] - gravity darkening exponent

        supersample_factor (int): number of points subdividing exposure
        exp_time (float): Exposure time (in same units as `time`)
        """

        pass
