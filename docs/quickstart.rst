.. quickstart

Quickstart
==========
Below is an example of basic ``evilmc`` usage to calculate a model light curve
including just the ellipsoidal variations and beaming signals. See
:class:`~evilmc.evparams` for detailed description of system parameters.

::

    import numpy as np
    import matplotlib.pylab as plt

    from evilmc import evparams, evmodel
    time = np.linspace(0, 1., 100)

    ep = evparams(per=1., a=4.15, T0=0.5, p=1./12.85, 
        limb_dark='quadratic', u=[0.314709, 0.312125], beta=0.07, 
        q=1.10e-3, Kz=0., Ts=6350., Ws=[0.,0.,0.1])

    em = evmodel(time, ep, 
        supersample_factor=5, exp_time=np.max(time)/time.shape)

    signal = em.evilmc_signal(num_grid=31)
    plt.plot(time, signal, ls='', marker='.')
    plt.show()

