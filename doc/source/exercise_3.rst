Exercise 3
----------

For the third exercise, the simple bi-directional scattering was replaced
by Rayleigh scattering and a more realistic, reflecting surface was added
at the bottom of the domain.

Similar as for the previous exercises, the simulations can be run by building
the :code:`exercise_3` target using :code:`make`:

.. code:: bash

 $ make exercise_3

Part 1
======

The first part of the exercise consisted of the simulation of a homogeneous
medium with the following properties

- Optical-depth: :math:`\tau` = 1
- Single scattering albedo: :math:`a` = 0.5
- Rayleight-scattering phase function
- Lower horizontal surface with an albedo :math:`A = 0.8` and Lambertian reflection
- Periodic vertical boundaries.

The simulation was performed in a three-dimensional domain with a horizontal
extent of :math:`10` km, a vertical extent of :math:`100` km and a depth of
:math:`1` km. The depth dimension is not relevant here since the scattering
plane for all scattering processes was restricted to the :math:`xy`-plane, which
is why the simulation was effectively performed as a 2D simulation. For the
described medium the corresponding optical properties are:

- absorption coefficient :math:`\sigma_a = 0.5 \cdot 10^{-5}\ \text{m}^{-1}` 
- scattering coefficient :math:`\sigma_s = 0.5 \cdot 10^{-5}\ \text{m}^{-1}`

Benchmark results
~~~~~~~~~~~~~~~~~

The total absorbed intensity for the different parts of the simulated domain
are given in the table below. They deviate from the reference values
but remain within the margin of uncertainty considering the large reported
differences between the "simple" and the "exact" integration schemes.

+---------------+----------+-----------------+
| Lower surface |  Medium  |  Upper boundary | 
+===============+==========+=================+
|        0.104  | 0.661    |          0.236  |
+---------------+----------+-----------------+

Part 2
======

Varying surface albedo
~~~~~~~~~~~~~~~~~~~~~~

The results for varying the albedo :math:`A` of the lower surface are are displayed in
the figure below. The plot displays the fractions of the total incoming intensity that
are absorbed by the lower surface, the medium and the leave the medium through the
upper boundary. Quite naturally, the intensity absorbed by the surface decreases as
its albedo increases. As the albedo increases, both the intensity of the medium and
the intensity leaving the domain through the upper boundary increase. This is expected
since, by conservation of energy, photons that are not absorbed by the surface must
by absorbed by the medium or leave the domain through the upper boundary.

.. figure:: ../../bin/results_3_b_1.png
   :alt: Results of exercise 3 b.1

   Absorbed intensity in different regions of the domain as the surface albedo is varied.


Varying single scattering albedo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The results for varying the single scattering albedo :math:`a` of the medium  are displayed
in the figure below. While the fraction of intensity absorbed by the medium decreases as the
the single scattering albedo decreases, the fraction of upwelling intensity increases
significantly. The intensity absorbed by the surface does not change significantly but
exhibits a slightly convex shape.

.. figure:: ../../bin/results_3_b_2.png
   :alt: Results of exercise 3 b.1

   Absorbed intensity in different regions of the domain as the the single scattering albedo
   is varied.

Varying the optical depth
~~~~~~~~~~~~~~~~~~~~~~~~~

The results for varying the optical depth :math:`\tau` of the medium are
displayed in the figure below. As could be expected, the fraction of intensity
absorbed by the medium increases as the optical depth is increased. At the same time
the intensity of the upwelling radiation is decreased. The intensity absorbed by the
surface decreases slightly.

.. figure:: ../../bin/results_3_b_3.png
   :alt: Results of exercise 3 b.3

   Absorbed intensity in different regions of the domain as the optical depth is varied.
