# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import simps
from scipy.misc import derivative

from astropy import constants as const

# https://stackoverflow.com/questions/779495/python-access-data-in-package-subdirectory
import pkg_resources

from PyAstronomy.modelSuite.XTran.forTrans import MandelAgolLC
from PyAstronomy import pyasl
import transit_utils

__all__= ['evmodel', 'evparams']

class evmodel(object):
    """Returns ellipsoidal variation of a slowly-rotating star induced by a 
    low-mass companion

    Args:
        time (numpy array): observational time (same units at period)
        params (:attr:`evparams`): object containing system parameters
        supersample_factor (int, optional): 
            number of points subdividing exposure
        exp_time (float, optional): Exposure time (in same units as `time`)
        response_function (str, optional): "Kepler" or "TESS";
            defaults to Kepler (and TESS isn't implemented yet)
    """

    def __init__(self, time, params,
            supersample_factor=1, exp_time=0, 
            response_function="Kepler"):
        """__init__ method for EVILMC
        """
        self.time = time
        self.params = params

        # Consider finite exposure time
        self.supersample_factor = supersample_factor
        self.exp_time = exp_time
        self.time_supersample = time
        if(self.exp_time > 1):
            self.time_supersample =\
                    transit_utils.supersample_time(time, supersample_factor,
                            exp_time)

        self.response_function = response_function

        # Calculate orbital phase
        phase_params = {"per": params.per, "T0": params.T0}
        self.phase = transit_utils.calc_phi(self.time_supersample, 
                phase_params)

    def evilmc_signal(self, num_grid=31):
        """Calculates the ellipsoidal variation and beaming effect curves

        Args:
            num_grid (int, optional): # of lat/long grid points on star

        Returns:
            numpy array: time-series ellipsoidal variation and beaming signals
        """
        # Make grid on stellar surface
        grid = _stellar_grid_geometry(self.params, num_grid)

        # Calculate orbital inclination in degrees
        inc = np.arccos(self.params.b/self.params.a)*180./np.pi

        # Calculate 3D orbital position of companion
        ke = pyasl.KeplerEllipse(self.params.a, self.params.per, 
                i=inc, tau=self.params.T0)
        rc = ke.xyzPos(self.time_supersample)

        # Calculate radial distance between companion and host
        nrm_rc = ke.radius(self.time_supersample)

        # Unit vector pointing toward planet
        rc_hat = rc/nrm_rc[:, None]

        # z-projection of orbital velocity, in fractions of speed of light
        vz = ke.xyzVel(self.time_supersample)[:, 2]

        # For each point in the orbit,
        for i in range(len(vz)):

            # Because the variation in stellar surface temperature with gravity
            # is so small, we approximate the variation in stellar radiation
            # using a first-order Taylor expansion. 
            # This approach also speeds up the calculation.
    
            # Calculate stellar radiation 
            # convolved with response function and Doppler shifts
            strad = _calc_stellar_brightness(
                    self.response_function, 
                    self.params.Ts,
                    vz[i])
            # Calculate temperature derivative of stellar radiation
            wrapped = lambda x:\
                    _calc_stellar_brightness(
                            self.response_function, 
                            x, 
                            vz[i])
            dx = self.params.Ts/1000.
            dstrad = derivative(wrapped, self.params.Ts, dx=dx)

            # psi is the angle between the companion's position vector
            # and the position vector of the center of the stellar grid element
            

        # Calculate the deformation for a very slightly tidally deformed 
        # and slowly rotating body with a Love number of 1

        # Calculate the small correction to the surface gravity vector 
        # for a very slightly tidally deformed and slowly rotating body

        return None

def _calc_stellar_brightness(which_response_function, Ts, vz):
    """Conolves the stellar radiation model 
    with the instrument response function
    at given temperature and for a given Doppler velocity

    Args:
        Ts (float): stellar temperature (K)
        vz (float): Doppler velocity in fractions of light speed
        which_response_function (str): "Kepler"

    """

    # MKS constants
    c = const.c.to('m/s').value
    h = const.h.to("J*s").value
    k_B = const.k_B.to("J/K").value

    if(which_response_function == "Kepler"):
        # https://keplergo.arc.nasa.gov/kepler_response_hires1.txt
        response_function_file =\
                pkg_resources.resource_filename('evilmc', 
                        'data/kepler_response_hires1.txt')
        # 2018 Jun 11 - This doesn't work!
        print(response_function_file)
    
        wavelength, resp =\
                np.genfromtxt(response_function_file,\
                comments="#", delimiter="\t", unpack=True)

        # Convert wavelengths from nanometers to meters
        wavelength *= 1e-9

        # Normalize
        resp /= simps(resp, wavelength)

        # Make into frequencies
        freq = c/wavelength[::-1]
    else:
        raise ValueError("which_response_function must be 'Kepler'!")

    # Using expression from Loeb & Gaudi (2003) ApJL 588, L117.
    freq0 = freq*(1. + vz)
    x0 = h*freq0/(k_B*Ts)
    # From Loeb & Gaudi, Eqn 3
    alpha0 = (np.exp(x0)*(3. - x0) - 3.)/(np.exp(x0) - 1.)

    F_nu0 = 2.*h*(freq0*freq0*freq0)/(c*c)/(np.exp(x0) - 1.)
    F_nu = F_nu0*(1. - (3. - alpha0)*vz)

    func = F_nu*resp

    return simps(func, freq)

class _stellar_grid_geometry(object):
    """Generates geometry for the stellar hemisphere facing the observer, 
    i.e. z > 0

    Args:
        params (:attr:`evparams`): object containing system parameters
    """

    def __init__(self, params, num_grid):

        # cos(theta) runs from 1 to 0 on the front face of the star,
        # and so the grid spacing dcos_theta is 1 over the grid number
        dcos_theta = 1./num_grid

        cos_theta = np.linspace(0.5*dcos_theta, 1. - 0.5*dcos_theta, num_grid)
        sin_theta = np.sqrt(1. - cos_theta**2.)

        # phi runs from 0 to 2 pi
        dphi = 2.*np.pi/(num_grid)
        phi = np.linspace(0.5*dphi, 2.*np.pi - 0.5*dphi, num_grid)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)

        # rhat is the normal vector to the stellar surface
        #
        # x/y/zhat are the x/y/z components of rhat and are defined on the
        # 2D grid of stellar surface grid points
        xhat = np.outer(cos_phi, sin_theta)
        yhat = np.outer(sin_phi, sin_theta)
        zhat = np.outer(np.ones_like(cos_phi), cos_theta)

        # nrm_Omega is the length of the stellar rotation vector
        Omega = params.Ws
        nrm_Omega = np.sqrt(
                Omega[0]*Omega[0] +\
                Omega[1]*Omega[1] +\
                Omega[2]*Omega[2]
                )

        # Omegahat is a unit vector pointing along Omega
        # If nrm_Omega is zero, then Omegahat is arbitrary
        if(nrm_Omega == 0.):
            Omegahat = np.array([0., 0., 1.0])
        else:
            Omegahat = Omega/nrm_Omega

        # cos_lambda is the cosine of the angle between rhat and the 
        # stellar rotation axis
        cos_lambda =\
                Omegahat[0]*xhat +\
                Omegahat[1]*yhat +\
                Omegahat[2]*zhat


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
        >>> ev = evparams(per=1., a=4.15, T0=0.5, p=1./12.85,
                limb_dark='quadratic', u=[0.314709, 0.312125], beta=0.07,
                q=1.10e-3, Kz=1e-6, Ts=6350., Ws=[0.,0.,1.])
        >>> # Print one example
        >>> print(ev.time)
    """

    # From https://stackoverflow.com/questions/8187082/how-can-you-set-class-attributes-from-variable-arguments-kwargs-in-python
    def __init__(self, **kwargs):

        # all those keys will be initialized as class attributes
        allowed_keys = set(['per', 'a', 'T0', 
            'p', 'limb_dark', 'u', 'beta', 'b', 'q', 'Kz', 'Ts', 'Ws'])
        # initialize all allowed keys to false
        self.__dict__.update((key, 0.) for key in allowed_keys)
        # and update the given keys by their given values
        self.__dict__.update((key, value) for key, value in kwargs.items()
                if key in allowed_keys)

