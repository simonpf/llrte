Exercise 1
----------

The first exercise consists of two simulations of radiative transfer
through purely absorbing media. The first one is restricted to a homogeneous
medium, while in the second one the medium is heterogeneous.

The simulations for the first exercise can be run by building the :code:`exercise_1`
target using :code:`make`:

.. code:: bash

 $ make exercise_1

Part 1
======

The first exercise was to simulate radiation propagating through a homogeneous,
purely absorbing medium.The simulated domain extends 10 km in the x-direction
and 1 km in the y- and z-direction. The source is a perfectly colimated beam
emitting photons in the x-direction. Absorption is set to a fixed value of
:math:`\alpha = 10^{-4}` to yield an optical depth of :math:`\tau = 1` over the simulated
domain.

Two experiments were conducted here. The first one investigates the impact of
the resolution of the grid on the results and the second on the impact of the
number of simulated photons.

Impact of grid resolution
~~~~~~~~~~~~~~~~~~~~~~~~~

Shown in the figure below are the fractions of photons that reach a given
distance from the source located at 0 km. Each set of results was obtained
by simulating :math:`10^4` photons. Colored markers show the results
for different grid resolutions. The black, dashed line shows the reference
result given by the absorption law. Since the domain is homogeneous, the
grid resolution does not have any significant impact on the results of the
simulation.

After all, since the domain is homogeneous it is not very surprising that
the grid resolution does not have a strong impact on the results.

.. figure:: ../../bin/results_1_a_1.png
   :alt: Results of exercise 1 a

   Monte carlo simulation results of a homogeneous purely absorbing medium for
   varying grid resolution.

Impact of photon count
~~~~~~~~~~~~~~~~~~~~~~

The impact of the number of simulated photons is shown in the figure below. Here
a grid resolution of 10 cells per kilometer was used and the number of simulated
photons was varied. For photon counts as low as 100 the results exhibit clearly
visible, random deviations from the reference curve. For photon counts higher
than that, however, they remain virtually indistinguishable.

.. figure:: ../../bin/results_1_a_2.png
   :alt: Results of exercise 1 a

   Monte carlo simulation results of a homogeneous purely absorbing medium for
   varying photon counts.


Part 2
======

For the second part, a heterogeneous domain has been simulated.


Implementation
~~~~~~~~~~~~~~

In the current MC implementation, the absorption model that determines the
absorption at each grid point is template parameter of the atmosphere to
simulate. The implementation from the first part can therefore easily be extended
to a heterogeneous domain by defining a custom absorption model class:

.. code:: c

 template <typename Float>
 // Absorption model for a heterogeneous domain consisting of two
 // parts with different absorption coefficient with a discrete
 // boundary along the x-axis.
 class HeterogeneousAbsorption {
 public:
     HeterogeneousAbsorption(Float abs_1,
                             Float abs_2,
                             Float boundary) :
         abs_1_(abs_1),
         abs_2_(abs_2),
         boundary_(boundary) {
         // Nothing to do.
     }
     // Called by MC simulation to get abs. coeff. at grid point.
     template <typename Grid, typename Position>
     Float get_absorption_coefficient(Grid /*grid*/, Position position) {
         if (position.x < boundary_) {
             return abs_1_;
         } else {
             return abs_2_;
         }
     }
 private:
     Float abs_1_, abs_2_, boundary_;
 };

Which is then used as absorption model in the simulated atmosphere:

.. code:: c

 auto absorption_model = HeterogeneousAbsorption<Float>(0.5e-4, 1.5e-4, 5e3);
 auto scattering_model = llrte::NoScattering<Float>();
 auto atmosphere = Atmosphere{grid, absorption_model, scattering_model};


Results
~~~~~~~

The configuration simulated here is an absorption coefficent :math:`\alpha = 0.5 \cdot 10^{-4}`
for :math:`x \leq 5` km and :math:`\alpha = 1.5 \cdot 10^{-4}` for the rest of the domain.
The results are shown for varying grid resolution and photon counts in the figures shown below.

Since the optical depth of the simulated domain is the same, the intensity at the far end of it
are the same as in the first simulation. Nonetheless, the effect of the heterogeneity of the
domain is clearly visible through the increased decrease of the intensity in the second half
of the domain. The qualitative results regarding grid resolution and photon count are the same
as for the first part.

.. figure:: ../../bin/results_1_b_1.png
   :alt: Results of exercise 1 b

   Monte carlo simulation results for a heterogeneous, purely absorbing medium for
   varying grid resolution.

.. figure:: ../../bin/results_1_b_2.png
   :alt: Results of exercise 1 b

   Monte carlo simulation results for a heterogeneous, purely absorbing medium for
   varying photon counts.
