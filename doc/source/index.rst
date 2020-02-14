.. llrte documentation master file, created by
   sphinx-quickstart on Fri Oct 18 07:28:42 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Monte Carlo methods for radiative transfer
==========================================

This project is an exploration of Monte Carlo techniques applied
to solve the radiative transfer equation (RTE).

Formulation
-----------

Homogeneous and purely absorbing medium
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Starting with the simplest case, a homogeneous, purely absorbing medium,
the RTE can be solved using the following Monte Carlo approach:

1. Create a photon at the source. Its optical depth :math:`\tau_\text{photon}` is zero.
2. Sample a random optical depth :math:`\tau` from an exponential distribution :math:`p(x) = \exp(-x)`
3. Integrate the absorption coefficient :math:`\alpha` along the photon path to obtain the photon's optical
   depth :math:`\tau_\text{photon}`.
4. When :math:`\tau_\text{photon} = \tau` the photon is absorbed. Continue with the next photon.

It is easy to see that the fraction of photons reaching a given distance from their emission
source satisfies the RTE for the spectral intensity. For this, let :math:`X` denote the path that
a photon travels prior to its absorption. Since :math:`\tau_\text{photon}` follows an exponential
distribution, its CDF is given by :math:`F(x) = 1 - \exp(-x)`. Noting that this is simply
the complement of the fraction of photons traveling at least a distance :math:`\tau` it follows:

.. math::

   P(\tau_\text{photon} > \tau) = \exp(- \tau)

Up to boundary conditions, this is equivalent to the extinction law, which is the solution of
the RTE for a homogeneous, purely absorbing medium.

Scattering
~~~~~~~~~~

The above algorithm can be easily extended to allow for multiple scattering. All that
is required for this is to 

1. Create a photon at the source. Its optical depth :math:`\tau_\text{photon}` is zero.
2. Sample a random optical depth :math`tau` from an exponential distribution :math:`p(x) = \exp(-x)`
3. Integrate the attenuation coefficient :math:`\mu` along the photon path to obtain the photon's optical
   depth :math:`\tau_\text{photon}`.
4. If :math:`\tau_\text{photon} \geq \tau`, sample another random number :math:`r`

   a) If :math:`r` is smaller or equal than the single scattering albedo, the photon is scattered. Sample
      a new direction from the phase function and continue with 2.
   b) If :math:`r` is larger than the single scattering albedo, the photon is absorbed. Continue
      with the next photon.

Alternatively, absorption can be treated in a non-random manner. For this a
photon energy is introduced, which is set to :math:`1.0` upon creation of the photon
and which is reduced at each step by a factor of :math:`\exp(-\int_{\Delta s} \alpha\ ds)`.
In this case the absorption coefficient is neglected in the
calculation of :math:`\tau_\text{photon}` and is just the path integral of the scattering
coefficient. The photon is then tracked until its energy decreases below a
suitably chosen threshold.

Implementation
--------------

:code:`llrte` is implemented in C++. So far it can perform simple non-polarized RT
calculations in 3 dimensions. It makes fairly heavy use of template meta
programming with the aim of keeping the implementation general at the same time
as being computationally efficient. Another reason is of course to scare
away FORTRAN programmers.

Exercises
---------

The solutions to the specific solution exercises can be found below.

.. toctree::
    :maxdepth: 1

    ./exercise_1
    ./exercise_2
    ./exercise_3
    ./exercise_4
    ./exercise_5
