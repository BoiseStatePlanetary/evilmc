# -*- coding: utf-8 -*-

import numpy as np

from PyAstronomy.modelSuite.XTran.forTrans import MandelAgolLC
from PyAstronomy import pyasl
import transit_utils

__all__= ['evmodel']

class evmodel(object):
    """Returns ellipsoidal variation of a slowly-rotating star induced by a 
    low-mass companion

    Args:
        time (numpy array): observational time (same units at period)
        params (:attr:`evparams`): object containing system parameters
        supersample_factor (int, optional): 
            number of points subdividing exposure
        exp_time (float, optional): Exposure time (in same units as `time`)
        num_grid (int, optional): # of lat/long grid points on star, 
            defaults to 31
    """

    def __init__(self, time, params, supersample_factor=1, exp_time=0, 
            num_grid=31):
        """__init__ method for EVILMC
        """
        self.time = time
        self.params = params

        self.supersample_factor = supersample_factor
        self.exp_time = exp_time

        self.time_supersample = time
        if(self.exp_time > 1):
            self.time_supersample =\
                    transit_utils.supersample_time(time, supersample_factor,
                            exp_time)

        # Calculate orbital phase
        phi_params = {"per": params.per, "T0": params.T0}
        self.phi = transit_utils.calc_phi(time, phi_params)

        # Calculate 3D orbital position of planet


class evparams(object):
    """System parameters for EVIL-MC calculation

    Args:
        per (float): orbital period (any units).
        a (float): semi-major axis (units of stellar radius)
        T0 (float):  mid-transit time (same units as period)
        p (float): planet's radius (units of stellar radius)
        limb_dark (str): Limb darkening model
            (choice of 'nonlinear' or 'quadratic')
        u (array_like): list of limb-darkening coefficients
        beta (float): gravity darkening exponent
        b (float): impact parameter (units of stellar radius)
        q (float): planet-star mass ratio
        Kz (float): stellar reflex velocity amplitude (in lightspeed)
        Ts (float): stellar temperature (in K)
        Ws (numpy array): stellar rotation vector in x, y, z

    Example:
        >>> import numpy as np
        >>> from evilmc import evparams
        >>> time = np.linspace(0, 1, 100)
        >>> ev = evparams(time=time, per=1., a=4.15, T0=0.5, p=1./12.85,
                limb_dark='quadratic', u=[0.314709, 0.312125], beta=0.07,
                q=1.10e-3, Kz=1e-6, Ts=6350., Ws=[0.,0.,1.])
        >>> # Print one example
        >>> print(ev.time)
    """

    # From https://stackoverflow.com/questions/8187082/how-can-you-set-class-attributes-from-variable-arguments-kwargs-in-python
    def __init__(self, **kwargs):

        # all those keys will be initialized as class attributes
        allowed_keys = set(['time','per','a', 'T0', 'p', 'limb_dark', 
            'u', 'beta', 'b', 'q', 'Kz', 'Ts', 'Ws'])
        # initialize all allowed keys to false
        self.__dict__.update((key, None) for key in allowed_keys)
        # and update the given keys by their given values
        self.__dict__.update((key, value) for key, value in kwargs.items()
                if key in allowed_keys)

