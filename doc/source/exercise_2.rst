Exercise 2
----------

For the second exercise, a domain was considered that would not only
absorb radiation, but also scatter it in forward or backward direction.
Again, the exercise consisted of two parts. In the first one a homogeneous
domain was to be simulated and in the second on a heterogeneous domain.

Similar as for the first exercise, the simulations can be run by building the :code:`exercise_2`
target using :code:`make`:

.. code:: bash

 $ make exercise_2

Part 1
======

The first part of the exercise consisted of the simulation of a homogeneous
medium with the following properties

- Optical-depth: :math:`\tau` = 1
- Single scattering albedo: :math:`a` = 0.8
- Forward-to-backscattering ratio :math:`VR = 0.8`

For a homogeneous slab of an extent of 10 km, this translates to an
absorption coefficient :math:`\sigma_a = 0.2 \cdot 10^{-4}\ \text{m}^{-1}` and
scattering coefficient :math:`\sigma_s = 0.8 \cdot 10^{-4}\ \text{m}^{-1}`.

Benchmark results
~~~~~~~~~~~~~~~~~

The following values were observed for the amounts of absorbed energy in the
domain versus the amount of energy leaving the domain through the boundaries:

+-----------------+----------+-----------------+
| Left boundary   | Medium   | Right boundary  |
+=================+==========+=================+
|        0.1172   | 0.1801   |          0.702  |
+-----------------+----------+-----------------+

Varying grid resolution
~~~~~~~~~~~~~~~~~~~~~~~

The results for varying grid sizes are displayed in the figure below. Compared
to the first exercise, the high single-scattering albedo and forward-to-backscattering
ratio effectively reduces amount of radiation the is absorbed by the domain.

A quite drastic effect of the grid resolution on the amount of scattered intensity
is visible in the results. This is not surprising since the number of samples per
grid cell decreases as the number of grid cells decreases. For the absorption, however,
this effect is not observed since photons continuously deposit energy along their path
thus increasing the number of samples per bin.

.. figure:: ../../bin/results_2_a_1.png
   :alt: Results of exercise 2 a

   Monte carlo simulation results of a homogeneous absorbing and scattering medium for
   varying grid resolution.

Varying photon counts
~~~~~~~~~~~~~~~~~~~~~

For variations in the number of simulated photons an effect on both the absorbed intensity
as well as the scattered intensity is observed. As expected, for low photon counts, the results
become more noisy.

.. figure:: ../../bin/results_2_a_2.png
   :alt: Results of exercise 2 a

   Monte carlo simulation results of a homogeneous absorbing and scattering medium for
   varying photon counts.

Varying optical depth
~~~~~~~~~~~~~~~~~~~~~

Quite naturally, the results become more interesting as the optical properties of the medium
are varied. The figure below shows the how the radiative transfer through the medium is
modified as its optical depth is varied. These results were obtained with 100 grid cells and
:math:`10^4` sampled photons. Quite intuitively, the absorption as well as the scattering
increases as the optical depth of the medium is increased.

.. figure:: ../../bin/results_2_a_3.png
   :alt: Results of exercise 2 a

   Monte carlo simulation results of a homogeneous absorbing and scattering medium for
   varying the optical depth of the medium.

Varying forward-to-backscattering ratio
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The next figure displays the effect of varying the forward-to-backscattering
ratio between :math:`0.2` and :math:`0.8`. The effect of this is most clearly
visible in absorption curves show in panel (a). A decreased value of :math:`VR`
leads to increased absorption closer to the source and decreased absorption
further away from the source. For high values of :math:`VR` the effect is the
opposite. Overall a lower forward-to-backscattering ratio increases the amount
of radiation absorbed by the medium. This is a fairly intuitive result since
backscattering event increase the path length of photons traveling through the
medium.

.. figure:: ../../bin/results_2_a_4.png
   :alt: Results of exercise 2 a

   Monte carlo simulation results of a homogeneous absorbing and scattering medium for
   varying the forward-to-backscattering ratio.

Varying single-scattering albedo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finally, the effect of varying the single-scattering albedo is investigated. The results
are displayed in the figure below. A lowered single-scattering albedo increases the absorption
in the medium, while an increased single-scattering albedo decreases it. Also this result
is quite intuitive since the lower single scattering albedo corresponds to an increased
absorption cross section and considering the high foward-to-backscattering ratio of the medium.

.. figure:: ../../bin/results_2_a_5.png
   :alt: Results of exercise 2 a

   Monte carlo simulation results of a homogeneous absorbing and scattering medium for
   varying the single-scattering albedo.

Part 2
======

We briefly show results for a heterogeneous medium. The simulation uses the same
geometry as above with an optical depth :math:`\tau = 1`. The heterogeneity is
introduced at :math:`x = 5` km, where the single scattering changes from
:math:`VR = 0.8` to :math:`VR = 0.2`. The results are shown in the figure below.
The effect of the increased absorption caused by the jump in forward to
back-scattering ratio is clearly visible in the results. In particular, the
non-local effect of the scattering properties of the medium become apparent,
since the absorption is increased even in the first half of the domain.
Regarding the frequency of scattering events, it is noticeable that the number
of zero scattering event is unchanged by the heterogeneity, which is expected as
the scattering cross section of the material is unchanged.

.. figure:: ../../bin/results_2_b.png
   :alt: Results of exercise 2 a

   Monte carlo simulation results of a heterogeneous absorbing and scattering medium.

Part 3
======

For the final part of the exercise, a reflecting surface was added at the end of
the simulation domain. The absorbed intensity in the different parts of the
domain are shown in the table. They agree well with the benchmark values.

+-----------------+----------+-----------------+
| Left boundary   | Medium   | Right boundary  |
+=================+==========+=================+
|        0.299    | 0.217    |          0.484  |
+-----------------+----------+-----------------+

The plot of the simulation results displays an increase of the absorbed as well
as the scattered intensity over the whole domain. This is easily understood
considering that photons that are reflected by the surface travel a longer
distance through the medium and are thus more likely to be absorbed or
scattered.


.. figure:: ../../bin/results_2_c.png
   :alt: Results of exercise 2 a

   Monte carlo simulation of a homogeneous medium with reflection.
