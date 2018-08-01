# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import simps
from scipy.misc import derivative
from scipy.integrate import simps

from astropy import constants as const
from astropy import units as u

from PyAstronomy.modelSuite.XTran.forTrans import MandelAgolLC

import pkg_resources

# MKS constants
c = const.c.to('m/s').value
h = const.h.to("J*s").value
k_B = const.k_B.to("J/K").value

__all__= ['evmodel', 'evparams', 'convert_Kz']

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

        # Calculate orbital inclination in degrees
        self.params.inc = np.arccos(self.params.b/self.params.a)*180./np.pi

        # Consider finite exposure time
        self.supersample_factor = supersample_factor
        self.exp_time = exp_time
        self.time_supersample = time
        if(self.supersample_factor > 1):
            self.time_supersample =\
                    _supersample_time(time, 
                            self.supersample_factor, self.exp_time)

        self.which_response_function = which_response_function
        self.response_function =\
                _retreive_response_function(self.which_response_function)

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

    def all_signals(self, num_grid=31):
        """Returns all signals, transits and eclipses included

        Args:
            num_grid (int, optional): # of lat/long grid points on star

        Returns:
            numpy array: time-series 

        Example:
            >>> import numpy as np
            >>> import matplotlib.pylab as plt
            >>> from evilmc import evparams, evmodel
            >>> time = np.linspace(0, 1., 100)
            >>> ep = evparams(per=1., a=4.15, T0=0.5, p=1./12.85,
            >>>     limb_dark='quadratic', u=[0.314709, 0.312125], beta=0.07,
            >>>     q=1.10e-3, Kz=0., Ts=6350., Ws=[0.,0.,0.1],
            >>>     F0=30e-6, Aplanet=30e-6, phase_shift=0.)
            >>> em = evmodel(time, ep,
            >>>     supersample_factor=5, exp_time=np.max(time)/time.shape)
            >>> signal = em.all_signals(num_grid=31)
            >>> plt.plot(time, signal, ls='', marker='.')
            >>> plt.show()
        """

        # Because I want to zero-out the eclipse portion of the signal
        #   even if the eclipse isn't included in the original request,
        #   I'm tacking on one time point, right in the middle of the 
        #   eclipse.
        TE = _calc_eclipse_time(self.params)
        self.time_supersample = np.append(self.time_supersample, TE)

        # Calculate orbital phase
        phase_supersample = _calc_phi(self.time_supersample, self.params)
        
        transit = self._transit() - 1.

        eclipse_depth = self.params.F0 + self.params.Aplanet
        eclipse = self._eclipse(eclipse_depth)

        E = self._calc_evilmc_signal(num_grid=num_grid)

        F0 = self.params.F0
        Aplanet = self.params.Aplanet
        phase_shift = self.params.phase_shift
        R = _reflected_emitted_curve(phase_supersample,\
                F0, Aplanet, phase_shift)

        ret = transit + E + R*eclipse

        # Now remove the extra time point at the end of the array
        shift_point = ret[-1]
        ret -= shift_point
        ret = ret[:-1]
        self.time_supersample = self.time_supersample[:-1]

        # Downsample if necessary
        if(self.supersample_factor > 1):
            ret = np.mean(ret.reshape(-1, self.supersample_factor), axis=1)

        return ret

    def evilmc_signal(self, num_grid=31):
        """Calculates the ellipsoidal variation and beaming effect curves

        Args:
            num_grid (int, optional): # of lat/long grid points on star

        Returns:
            numpy array: time-series ellipsoidal variation and beaming signals

        Example:
            >>> import numpy as np
            >>> import matplotlib.pylab as plt
            >>> from evilmc import evparams, evmodel
            >>> time = np.linspace(0, 1., 100)
            >>> ep = evparams(per=1., a=4.15, T0=0.5, p=1./12.85,
            >>>     limb_dark='quadratic', u=[0.314709, 0.312125], beta=0.07,
            >>>     q=1.10e-3, Kz=0., Ts=6350., Ws=[0.,0.,0.1])
            >>> em = evmodel(time, ep, supersample_factor=5,
            >>>     exp_time=np.max(time)/time.shape)
            >>> signal = em.evilmc_signal(num_grid=31)
            >>> plt.plot(time, signal, ls='', marker='.')
            >>> plt.show()
        """

        # Because I want to zero-out the eclipse portion of the signal
        #   even if the eclipse isn't included in the original request,
        #   I'm tacking on one time point, right in the middle of the 
        #   eclipse.
        TE = _calc_eclipse_time(self.params)
        self.time_supersample = np.append(self.time_supersample, TE)

        ret = self._calc_evilmc_signal(num_grid)

        # Now remove the extra time point at the end of the array
        shift_point = ret[-1]
        ret -= shift_point
        ret = ret[:-1]
        self.time_supersample = self.time_supersample[:-1]

        # Downsample if necessary
        if(self.supersample_factor > 1):
            ret = np.mean(ret.reshape(-1, self.supersample_factor),\
                    axis=1)

        return ret

    def _calc_evilmc_signal(self, num_grid):
        """Calculates the ellipsoidal variation and beaming effect curves

        Returns:
            numpy array: time-series ellipsoidal variation and beaming signals
        """
        # Make grid on stellar surface
        grid = _stellar_grid_geometry(self.params, self.Omegahat, num_grid)

        # Calculate 3D orbital position of companion
        rc = _calc_orbital_position(self.time_supersample, self.params)

        # Calculate radial distance between companion and host
        nrm_rc = self.params.a

        # Unit vector pointing toward planet
        rc_hat = rc/nrm_rc

        # z-projection of orbital velocity, in fractions of speed of light
        vz = self.params.Kz*_calc_projected_velocity(self.time_supersample, 
                                                        self.params)

        #integrated disk brightness
        stellar_disk = np.zeros_like(vz)

        # Because the variation in stellar surface temperature with gravity
        # is so small, we approximate the variation in stellar radiation
        # using a first-order Taylor expansion. 
        # This approach also speeds up the calculation.
    
        # Calculate stellar radiation 
        # convolved with response function and Doppler shifts
        strad = _calc_stellar_brightness(self.params.Ts, vz, 
                self.response_function)
        # Calculate temperature derivative of stellar radiation
        #
        # Since we're calculating the *temperature* derivative, the value of
        #   vz doesn't matter, so to avoid calling the function several times
        #   I'll just take the first value in vz.
        wrapped = lambda x:\
                _calc_stellar_brightness(x, vz[0], self.response_function)
        dx = self.params.Ts/1000.
        dstrad_dtemp = derivative(wrapped, self.params.Ts, dx=dx)

        cos_psi =   rc_hat[:,0]*grid.xhat[:,:,None] +\
                    rc_hat[:,1]*grid.yhat[:,:,None] +\
                    rc_hat[:,2]*grid.zhat[:,:,None] 

        # Calculate the deformation for a very slightly tidally deformed 
        # and slowly rotating body with a Love number of 1
        del_R = _calc_del_R(self.params.q, nrm_rc, cos_psi, self.nrm_Omega, 
                grid.cos_lambda)

        # Calculate the small correction to the surface gravity vector 
        # for a very slightly tidally deformed and slowly rotating body
        del_gam_vec_x = _del_gam_vec(del_R, grid.xhat, self.params.q, 
                nrm_rc, rc_hat[:, 0], cos_psi, self.nrm_Omega, 
                self.Omegahat[0], grid.cos_lambda)
        del_gam_vec_y = _del_gam_vec(del_R, grid.yhat, self.params.q, 
                nrm_rc, rc_hat[:, 1], cos_psi, self.nrm_Omega, 
                self.Omegahat[1], grid.cos_lambda)
        del_gam_vec_z = _del_gam_vec(del_R, grid.zhat, self.params.q, 
                nrm_rc, rc_hat[:, 2], cos_psi, self.nrm_Omega, 
                self.Omegahat[2], grid.cos_lambda)

        # x/y/z components of the local graviational acceleration
        gz = -grid.zhat[:,:,None] + del_gam_vec_z

        # dot product between rhat and the components of the 
        # gravity-correction vector
        rhat_dot_dgam = grid.xhat[:,:,None]*del_gam_vec_x +\
                grid.yhat[:,:,None]*del_gam_vec_y +\
                grid.zhat[:,:,None]*del_gam_vec_z

        # magnitude of modified local gravity vector
        nrm_g = 1. - rhat_dot_dgam

        # cos of angle between the line of sight and the gravity vector
        mu = abs(gz)/nrm_g

        dgam0 = _rhat_dot_del_gam0(self.params.q, nrm_rc, self.nrm_Omega)

        # small temperature correction
        dtemp = _del_temp(self.params.beta, rhat_dot_dgam, dgam0)*\
                self.params.Ts

        temp = self.params.Ts + dtemp
            
        # stellar radiation at temp
        strad_at_temp = np.ones_like(temp)*strad + dstrad_dtemp*dtemp

        #limb-darkened profile
        prof = _limb_darkened_profile(self.params.limb_dark, self.params.u,
                mu)

        # projected area of each grid element
        dareap = (1. + 2.*del_R)*mu*grid.dcos_theta*grid.dphi

        stellar_disk = np.sum(prof*strad_at_temp*dareap, axis=(0,1))
        stellar_disk = stellar_disk/np.nanmean(stellar_disk) - 1.

        return stellar_disk

    def _transit(self):
        """Uses PyAstronomy's quadratic limb-darkening routine to calculate
        the transit light curve

        Returns:
            numpy array: transit light curve
        """

        ma = MandelAgolLC(orbit="circular", ld="quad")

        # If quadratic limb-darkening
        if(self.params.limb_dark == 'quadratic'):
            ma["linLimb"] = self.params.u[0]
            ma["quadLimb"] = self.params.u[1]

        ma["per"] = self.params.per
        # Set using the impact parameter
        ma["i"] = self.params.inc
        ma["a"] = self.params.a
        ma["T0"] = self.params.T0
        ma["p"] = self.params.p

        return ma.evaluate(self.time_supersample)

    def _eclipse(self, eclipse_depth):
        """
        Uses PyAstronomy's transit light curve routine with uniform
        limb to calculate eclipse

        Args:
            eclipse_depth (float): eclipse depth
        """

        # Make eclipse_depth isn't zero!
        if(eclipse_depth != 0):
            ma = MandelAgolLC(orbit="circular", ld="quad")
            TE = _calc_eclipse_time(self.params)

            ma = MandelAgolLC(orbit="circular", ld="quad")

            # If quadratic limb-darkening
            if(self.params.limb_dark == 'quadratic'):
                ma["linLimb"] = 0.
                ma["quadLimb"] = 0.

            ma["per"] = self.params.per
            # Set using the impact parameter
            ma["i"] = self.params.inc
            ma["a"] = self.params.a
            ma["T0"] = TE
            ma["p"] = np.sqrt(eclipse_depth)

            eclipse = ma.evaluate(self.time_supersample)

            # Rescale eclipse
            eclipse = 1. - eclipse
            eclipse /= eclipse_depth
            eclipse = 1. - eclipse

        elif(eclipse_depth == 0.):
            eclipse = 0.

        return eclipse

def convert_Kz(Mp=None, Ms=None, q=None, a=0.1, inc=90., Kz=None):
    """Calculate and/or convert a radial velocity to fractions of the speed of
    light

    Args:
        Mp (float, optional): companion mass in units of Jupiter masses 
            (1.8981872e27 kg), defaults to None
        Ms (float, optional): host star mass in solar units 
            (1.9884754e30 kg), defaults to None
        q (float, optional): companion-host mass ratio, defaults to None
        a (float): semi-major axis in AU, defaults to 0.1
        inc (float): inclination angle in degrees, defaults to 90 degrees
        Kz (float, optional): radial velocity amplitude in m/s, defaults to
            None

    Returns:
        float: radial velocity as fraction of speed of light

    Example:
        >>> from evilmc import convert_Kz
        >>> Kz = 93. # in m/s - typical for short-period hot Jupiter
        >>> print(convert_Kz(Kz=Kz))
        
        >>> # Or if you have planetary and stellar mass
        >>> Mp = 1. # Jupiter masses
        >>> Ms = 1. # solar masses
        >>> print(convert_Kz(Mp=Mp, Ms=Ms))

        >>> # Or if you have planetary mass and the mass ratio
        >>> q = 1.e-3
        >>> Mp = 1.
        >>> print(convert_Kz(Mp=Mp, q=q))
    """
    
    MJup_to_MSun = (u.jupiterMass.to('kg'))/(u.solMass.to('kg'))

    def _calculate_RV(Mp, Ms, a, inc):

        # From Lovis & Fischer (2010), Eqn (13)
        RV0 = 28.4329/c #m/s to fractions of the speed of light
        return RV0*Mp*np.sin(inc*np.pi/180.)/\
                np.sqrt(Mp*MJup_to_MSun + Ms)/np.sqrt(a)

    if(Kz is not None):
        return Kz/c
    elif((Mp is not None) and (Ms is not None)):
        return _calculate_RV(Mp, Ms, a, inc)
    elif((Mp is not None) and (q is not None)):
        Ms = Mp/q*MJup_to_MSun
        return _calculate_RV(Mp, Ms, a, inc)

def _reflected_emitted_curve(phase, F0, Aplanet, phase_shift):
    """Returns sinusoidal reflection curve

    Args:
        phase (numpy array): orbital phase
        F0 (float): zero-point for phase curve
        Aplanet (float): phase curve amplitude
        phase_shift (float): phase curve shift
    """

    return F0 - Aplanet*np.cos(2.*np.pi*(phase - phase_shift))

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

    term0 = 2.*del_R*rhat[:,:,None]

    intermed_term = np.sqrt(a*a - 2.*a*cos_psi + 1.)
    term1 = q*(a*ahat -\
            rhat[:,:,None])/(intermed_term*intermed_term*intermed_term)

    term2 = nrm_Omega*nrm_Omega/(a*a*a)*(rhat-Omegahat*cos_lambda)[:,:,None]
    term3 = -q/(a*a)*ahat*np.ones_like(cos_lambda)[:,:,None]

    return term0 + term1 + term2 + term3

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
    term0 = q*(1./np.sqrt(r*r - 2.*r*cos_psi + 1.))
    term1 = -q*(1./np.sqrt(r*r + 1.) + cos_psi/(r*r))
    term2 = -nrm_Omega*nrm_Omega/(2.*r*r*r)*\
            (cos_lambda[:,:,None]*cos_lambda[:,:,None])

    return term0 + term1 + term2

def _retreive_response_function(which_response_function):
    """Retreives the instrument response function

    Args:
        which_response_function (str): "Kepler"

    Returns:
        dict of numpy arrays: "wavelength" (in meters) and "resp"
    """

    if(which_response_function == "Kepler"):
        # https://keplergo.arc.nasa.gov/kepler_response_hires1.txt
        path = 'data/kepler_response_hires1.txt'  # always use slash
        response_function_file = pkg_resources.resource_filename(__name__, path)

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

#   F_nu = F_nu0*(1. - (3. - alpha0)*vz)
    term0 = vz*(3. - alpha0)
    F_nu = F_nu0*(1. - term0)

    func = F_nu*resp[:,None]

    return np.trapz(func, freq, axis=0)

def _supersample_time(time, supersample_factor, exp_time):
    """Creates super-sampled time array

    Args:
        time (numpy array): times
        supersample_factor (int): number of points subdividing exposure
        exp_time (float): Exposure time (in same units as `time`)

    Returns:
        Returns the super-sampled time array
    """

    if supersample_factor > 1:
        time_offsets = np.linspace(-exp_time/2., exp_time/2.,
                supersample_factor)
        time_supersample = (time_offsets +\
                time.reshape(time.size, 1)).flatten()
    else:
        time_supersample = time

    return time_supersample

def _calc_phi(time, params):
    """Calculates orbital phase assuming zero eccentricity

    Args:
        time: observational time (same units at orbital period)
        params: dict of floats/numpy arrays, including
            params.per - orbital period (any units)
            params.T0 - mid-transit time (same units as period)

    Returns:
        orbital phase
    """

    T0 = params.T0
    per = params.per

    return ((time - T0) % per)/per

def _calc_eclipse_time(params):
    """Calculates mid-eclipse time - 
    I've included this function here in anticipation of using eccentric orbits 
    in the near future.
    """

    T0 = params.T0
    per = params.per

    return T0 + 0.5*per

def _calc_orbital_position(time, params):
    """Calculates x/y/z position of planet, 
    assuming zero eccentricity (as of 2018 Jul 30)

    Args:
        time: observational time (same units at orbital period)
        params (:attr:`evparams`): dict of floats/numpy arrays, including
            params.a - semi-major axis
            params.inc - orbital inclination
            params.T0 - mid-transit time
    """

    # mean motion
    n = 2.*np.pi/params.per

    # true anomaly
    f = 2.*np.pi*_calc_phi(time, params)

    # Ch. 2 of Murray & Dermott (1999), p. 51, Eqn 2.122
    # I'm assuming circular orbits here, 
    # so pericenter longitude assumed pi/2 and
    # ascending node longitude assumed 0 degrees
    xc = params.a*(np.cos(f + np.pi/2.))
    yc = params.a*(np.sin(f + np.pi/2.)*np.cos(params.inc))
    zc = params.a*(np.sin(f + np.pi/2.)*np.sin(params.inc))

    return np.array([xc, yc, zc]).transpose()

def _calc_projected_velocity(time, params):
    """Calculates the z-velocity of the host star

    Args:
        time: observational time (same units at orbital period)
        params (:attr:`evparams`): dict of floats/numpy arrays, including
            params.a - semi-major axis
            params.inc - orbital inclination
            params.T0 - mid-transit time
    """

    # mean motion
    n = 2.*np.pi/params.per

    # true anomaly
    f = 2.*np.pi*_calc_phi(time, params)

    # Remember that the amplitude of the velocity is a free parameter
    return np.sin(params.inc)*(np.cos(f + np.pi/2.))

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
        >>> ep = evparams(per=1., a=4.15, T0=0.5, p=1./12.85,
        >>>     limb_dark='quadratic', u=[0.314709, 0.312125], beta=0.07,
        >>>     q=1.10e-3, Kz=1e-6, Ts=6350., Ws=[0.,0.,0.1])
        >>> # Print one example
        >>> print(ev.per)
    """

    # From https://stackoverflow.com/questions/8187082/how-can-you-set-class-attributes-from-variable-arguments-kwargs-in-python
    def __init__(self, **kwargs):

        # all those keys will be initialized as class attributes
        allowed_keys = set(['per', 'a', 'T0', 
            'p', 'limb_dark', 'u', 'beta', 'b', 'q', 'Kz', 'Ts', 'Ws',
            'F0', 'Aplanet', 'phase_shift'])
        # initialize all allowed keys to false
        self.__dict__.update((key, 0.) for key in allowed_keys)
        # and update the given keys by their given values
        self.__dict__.update((key, value) for key, value in kwargs.items()
                if key in allowed_keys)

