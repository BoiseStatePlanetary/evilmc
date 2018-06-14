# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import simps
from scipy.misc import derivative

import matplotlib.pyplot as plt

from astropy import constants as const

from PyAstronomy.modelSuite.XTran.forTrans import MandelAgolLC
from PyAstronomy import pyasl
import transit_utils

# MKS constants
c = const.c.to('m/s').value
h = const.h.to("J*s").value
k_B = const.k_B.to("J/K").value

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
            which_response_function="Kepler"):
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

        self.which_response_function = which_response_function
        self.response_function =\
                _retreive_response_function(self.which_response_function)

        # Calculate orbital phase
        phase_params = {"per": params.per, "T0": params.T0}
        self.phase = transit_utils.calc_phi(self.time_supersample, 
                phase_params)

        # nrm_Omega is the length of the stellar rotation vector
        Omega = params.Ws
        self.nrm_Omega = np.sqrt(
                Omega[0]*Omega[0] +\
                Omega[1]*Omega[1] +\
                Omega[2]*Omega[2]
                )

        # Omegahat is a unit vector pointing along Omega
        # If nrm_Omega is zero, then Omegahat is arbitrary
        if(self.nrm_Omega == 0.):
            self.Omegahat = np.array([0., 0., 1.0])
        else:
            self.Omegahat = params.Ws/self.nrm_Omega


    def evilmc_signal(self, num_grid=31):
        """Calculates the ellipsoidal variation and beaming effect curves

        Args:
            num_grid (int, optional): # of lat/long grid points on star

        Returns:
            numpy array: time-series ellipsoidal variation and beaming signals
        """
        # Make grid on stellar surface
        grid = _stellar_grid_geometry(self.params, self.Omegahat, num_grid)

        # Calculate orbital inclination in degrees
        inc = np.arccos(self.params.b/self.params.a)*180./np.pi

        # Calculate 3D orbital position of companion
        ke = pyasl.KeplerEllipse(self.params.a, self.params.per, 
                i=inc, tau=self.params.T0, w=90.)
        rc = ke.xyzPos(self.time_supersample)

        # Calculate radial distance between companion and host
        nrm_rc = ke.radius(self.time_supersample)

        # Unit vector pointing toward planet
        rc_hat = rc/nrm_rc[:, None]

        # z-projection of orbital velocity, in fractions of speed of light
        vz = ke.xyzVel(self.time_supersample)[:, 2]/c

        # normalize and shift vz so Kz is a free parameter
        vz = vz/np.nanmax(vz)*self.params.Kz

        #integrated disk brightness
        disk = np.zeros_like(vz)

        # Calculate stellar radiation 
        # convolved with response function and Doppler shifts
        strad = _calc_stellar_brightness(self.params.Ts, vz, 
                self.response_function)
        # Calculate temperature derivative of stellar radiation
        wrapped = lambda x:\
                _calc_stellar_brightness(x, vz, self.response_function)
        dx = self.params.Ts/1000.
        dstrad_dtemp = derivative(wrapped, self.params.Ts, dx=dx)

        # For each point in the orbit,
        for i in range(len(vz)):

            # Because the variation in stellar surface temperature with gravity
            # is so small, we approximate the variation in stellar radiation
            # using a first-order Taylor expansion. 
            # This approach also speeds up the calculation.
    
            # psi is the angle between the companion's position vector
            # and the position vector of the center of the stellar grid element
            cos_psi =\
                    rc_hat[i, 0]*grid.xhat +\
                    rc_hat[i, 1]*grid.yhat +\
                    rc_hat[i, 2]*grid.zhat

            # Calculate the deformation for a very slightly tidally deformed 
            # and slowly rotating body with a Love number of 1
            del_R = _calc_del_R(self.params.q, nrm_rc[i], cos_psi, 
                    self.nrm_Omega, grid.cos_lambda)

            # Calculate the small correction to the surface gravity vector 
            # for a very slightly tidally deformed and slowly rotating body
            del_gam_vec_x = _del_gam_vec(del_R, grid.xhat, self.params.q, 
                    nrm_rc[i], rc_hat[i, 0], cos_psi, self.nrm_Omega, 
                    self.Omegahat[0], grid.cos_lambda)
            del_gam_vec_y = _del_gam_vec(del_R, grid.yhat, self.params.q, 
                    nrm_rc[i], rc_hat[i, 1], cos_psi, self.nrm_Omega, 
                    self.Omegahat[1], grid.cos_lambda)
            del_gam_vec_z = _del_gam_vec(del_R, grid.zhat, self.params.q, 
                    nrm_rc[i], rc_hat[i, 2], cos_psi, self.nrm_Omega, 
                    self.Omegahat[2], grid.cos_lambda)

            # x/y/z components of the local graviational acceleration
            gx = -grid.xhat + del_gam_vec_x
            gy = -grid.yhat + del_gam_vec_y
            gz = -grid.zhat + del_gam_vec_z

            # dot product between rhat and the components of the 
            # gravity-correction vector
            rhat_dot_dgam = grid.xhat*del_gam_vec_x +\
                    grid.yhat*del_gam_vec_y +\
                    grid.zhat*del_gam_vec_z
            # magnitude of modified local gravity vector
            nrm_g = 1. - rhat_dot_dgam

            # cos of angle between the line of sight and the gravity vector
            mu = abs(gz)/nrm_g

            dgam0 = _rhat_dot_del_gam0(self.params.q, nrm_rc[i], 
                    self.nrm_Omega)

            # small temperature correction
            dtemp = _del_temp(self.params.beta, rhat_dot_dgam, dgam0)*\
                    self.params.Ts

            temp = self.params.Ts + dtemp
            
            # stellar radiation at temp
            strad_at_temp = np.ones_like(temp)*strad[i] + dstrad_dtemp[i]*dtemp

            #limb-darkened profile
            prof = _limb_darkened_profile(self.params.limb_dark, self.params.u,
                    mu)

            # projected area of each grid element
            dareap = (1. + 2.*del_R)*mu*grid.dcos_theta*grid.dphi

            disk[i] = np.sum(prof*strad_at_temp*dareap)

        # normalize to minimum point
        disk /= np.nanmin(disk)

        return disk

def _limb_darkened_profile(limb_dark_law, LDCs, mu):
    """Returns limb-darkened flux

    Args:
        limb_dark_law (str): 'quadratic'
        LDCs (numpy array): limb-darkening coefficients
            Must have two elements for 'quadratic'
        mu (numpy array): cos of the angle between the line-of-sight and the
            normal to a grid element of the host's surface

    Returns:
        numpy array: normalized flux profile

    """

    if((limb_dark_law == 'quadratic') and (len(LDCs) == 2)):
        return 1. - LDCs[0]*(1. - mu) - LDCs[1]*(1. - mu)*(1. - mu)
    else:
        raise ValueError("Only quadratic limb-darkening law allowed, "+
                "which requires two coefficients!")

def _del_temp(beta, rhat_dot_dgam, dgam0):
    """Returns the small correction in surface temperature

    See Eqn (10) from Jackson+ (2012) ApJ.

    Args:
        beta (float): gravity-darkening exponent, probably 0.07 or 0.25
        rhat_dot_dgam (numpy array): dot product between the unit location
            vector for each grid element of the host's surface
        dgam0 (float): small correction to gravity vector at pole of host

    Returns:
        numpy array: small temperature correction at each point on the
            surface of the host

    """

    return beta*(dgam0 - rhat_dot_dgam)

def _rhat_dot_del_gam0(q, a, nrm_Omega):
    """Returns the small correction to the magnitude of the surface gravity
    vector at the 'pole' (see Jackson et al.) for a very slightly tidally 
    deformed and slowly rotating body

    See Eqn (9) from Jackson+ (2012) ApJ.

    Args: 
        q (float): companion-host mass ratio
        a (float): radial distance between the companion and host
        nrm_Omega (float): magnitude of the host's rotation vector

    Returns:
        float: correction to the polar gravity vector
    """

    term1 = np.sqrt(a*a + 1.)

    return -q/(term1*term1*term1) + nrm_Omega*nrm_Omega/(a*a*a)


def _del_gam_vec(del_R, rhat, q, a, ahat, cos_psi, nrm_Omega, Omegahat,
        cos_lambda):
    """Returns the small correction to the surface gravity vector for a very
    slightly tidally deformed and slowly rotating body
    
    See Eqn (8) from Jackson+ (2012) ApJ.

    Args: 
        del_R (numpy array): radial distance between the host's center and its
            surface for each element of the surface grid
        rhat (numpy array): unit vector pointing at the center of each grid
            element of the host's surface; probably x/y/zhat
        q (float): companion-host mass ratio
        a (float): radial distance between the companion and host
        ahat (float): x/y/z component of the normalized location vector for 
            the companion
        cos_psi (numpy array): cosine of the angle between the companion's
            mass location vector and the location vector for a point on the
            surface of the host
        nrm_Omega (float): magnitude of the host's rotation vector
        Omegahat (float): x/y/z component of the normalized rotation vector 
            for the host
        cos_lambda (numpy array): cosine of the angle between the
            companion's rotation vector and the location vector for a point 
            on the surface of the host
            

    Returns:
        numpy array: tiny deviation in gravitational vector components at 
            the center of each grid element for the host's surface
    """

    term1 = np.sqrt(a*a - 2.*a*cos_psi + 1.)

    return 2.*del_R*rhat +\
            q*(a*ahat - rhat)/(term1*term1*term1) +\
            nrm_Omega*nrm_Omega/(a*a*a)*(rhat-Omegahat*cos_lambda) -\
            q/(a*a)*ahat

def _calc_del_R(q, r, cos_psi, nrm_Omega, cos_lambda):
    """Returns the deformation for a very slightly tidally 
    deformed and slowly rotating body with a Love number of 1

    Args:
        q (float): companion-host mass ratio
        r (float): radial distance between the companion and host
        cos_psi (numpy array): cosine of the angle between the companion's
            mass location vector and the location vector for a point on the
            surface of the host
        nrm_Omega (float): magnitude of the host's rotation vector
        cos_lambda (numpy array): cosine of the angle between the
            companion's rotation vector and the location vector for a point 
            on the surface of the host

    Returns:
        numpy array: radial distance between the host's center and its 
            surface for each element of the surface grid
    """
    return q*(1./np.sqrt(r*r - 2.*r*cos_psi + 1.) -\
            1./np.sqrt(r*r + 1.) - cos_psi/(r*r)) -\
            nrm_Omega*nrm_Omega/(2.*r*r*r)*(cos_lambda*cos_lambda)

def _retreive_response_function(which_response_function):
    """Retreives the instrument response function

    Args:
        which_response_function (str): "Kepler"

    Returns:
        dict of numpy arrays: "wavelength" (in meters) and "resp"
    """

    if(which_response_function == "Kepler"):
        # https://keplergo.arc.nasa.gov/kepler_response_hires1.txt
        response_function_file = "data/kepler_response_hires1.txt"

        wavelength, resp =\
                np.genfromtxt(response_function_file,\
                comments="#", delimiter="\t", unpack=True)

        # Convert wavelengths from nanometers to meters
        wavelength *= 1e-9

        # Normalize
        resp /= simps(resp, wavelength)

    else:
        raise ValueError("which_response_function must be 'Kepler'!")

    return {"wavelength": wavelength, "resp": resp}



def _calc_stellar_brightness(Ts, vz, response_function):
    """Convolves the stellar radiation model 
    with the instrument response function
    at given temperature and for a given Doppler velocity

    Args:
        Ts (float): stellar temperature (K)
        vz (float): Doppler velocity in fractions of light speed
        response_function (dict of numpy arrays): "wavelength" and "resp"

    Returns:
        float: integrated host disk brightness in MKS; 
            exact value is unimportant since the signal time-series
            is normalized
    """

    wavelength = response_function['wavelength']
    resp = response_function['resp']

    # Make into frequencies
    freq = c/wavelength[::-1]

    # Using expression from Loeb & Gaudi (2003) ApJL 588, L117.
    freq0 = np.outer(freq, (1. + vz))
    x0 = h*freq0/(k_B*Ts)
    # From Loeb & Gaudi, Eqn 3
    alpha0 = (np.exp(x0)*(3. - x0) - 3.)/(np.exp(x0) - 1.)

    F_nu0 = 2.*h*(freq0*freq0*freq0)/(c*c)/(np.exp(x0) - 1.)
    F_nu = F_nu0*(1. - (3. - alpha0)*vz)

    func = (F_nu.transpose()*resp).transpose()

    return np.trapz(func, freq, axis=0)

class _stellar_grid_geometry(object):
    """Generates geometry for the stellar hemisphere facing the observer, 
    i.e. z > 0

    Args:
        params (:attr:`evparams`): object containing system parameters
        Omegahat (numpy array): x, y, and z components of the normalized 
            rotation vector for the host
    """

    def __init__(self, params, Omegahat, num_grid):

        # cos(theta) runs from 1 to 0 on the front face of the star,
        # and so the grid spacing dcos_theta is 1 over the grid number
        self.dcos_theta = 1./num_grid

        cos_theta = np.linspace(0.5*self.dcos_theta,
                1. - 0.5*self.dcos_theta, num_grid)
        sin_theta = np.sqrt(1. - cos_theta**2.)

        # phi runs from 0 to 2 pi
        self.dphi = 2.*np.pi/(num_grid)
        phi = np.linspace(0.5*self.dphi, 2.*np.pi - 0.5*self.dphi, num_grid)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)

        # rhat is the normal vector to the stellar surface
        #
        # x/y/zhat are the x/y/z components of rhat and are defined on the
        # 2D grid of stellar surface grid points
        self.xhat = np.outer(cos_phi, sin_theta)
        self.yhat = np.outer(sin_phi, sin_theta)
        self.zhat = np.outer(np.ones_like(cos_phi), cos_theta)

        # cos_lambda is the cosine of the angle between rhat and the 
        # stellar rotation axis
        self.cos_lambda =\
                Omegahat[0]*self.xhat +\
                Omegahat[1]*self.yhat +\
                Omegahat[2]*self.zhat


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
        Ws (numpy array): stellar rotation vector in x, y, z,
            in units of the orbital mean motion

    Example:
        >>> import numpy as np
        >>> from evilmc import evparams
        >>> ev = evparams(per=1., a=4.15, T0=0.5, p=1./12.85,
                limb_dark='quadratic', u=[0.314709, 0.312125], beta=0.07,
                q=1.10e-3, Kz=1e-6, Ts=6350., Ws=[0.,0.,0.118])
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

