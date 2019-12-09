Exercise 4
----------

For the fourth exercise, a medium with a Henyey-Greenstein phase function was
considered. The simulations can be run by building the :code:`exercise_4` target
using :code:`make`:

.. code:: bash

 $ make exercise_4

Part 1
======

The first part of the exercise consisted of the simulation of a homogeneous
medium with the following properties

- Optical-depth: :math:`\tau` = 5
- Single scattering albedo: :math:`a` = 0.8
- Henyey-Greenstein phase function with :math:`g` = 0.9
- Lower horizontal surface with an albedo :math:`A = 0.7` and Lambertian reflection
- Periodic vertical boundaries.

The simulation was performed in a three-dimensional domain with a horizontal
extent of :math:`4000` m, a vertical extent of :math:`100` m and a depth of
:math:`1` m. The depth dimension is not relevant here since the scattering
plane for all scattering processes was restricted to the :math:`xy`-plane, so
that the simulation was effectively performed as a 2D simulation. 
The corresponding optical properties for the described medium are:

- absorption coefficient :math:`\sigma_a = 1 \cdot 10^{-2}\ \text{m}^{-1}` 
- scattering coefficient :math:`\sigma_s = 4 \cdot 10^{-2}\ \text{m}^{-1}`

The Henyey-Greenstein phase function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Henyey-Greenstein phase function is given by the following expression

.. math::

   p(\theta) = \frac{1 - g^2}{(1 + g^2 + 2\ g\ \cos( \theta))^{\frac{3}{2}}}.

It is parametrized by the asymmetry parameter :math:`g`, which quantifies the
asymmetry of the scattering process with respect to the forward and backward
directions. The figure below displays the phase function for different values of
the asymmetry parameter :math:`g`. For negative :math:`g`, the curve has a
strong backward-scattering peak. As :math:`g` increases the backward scattering
peak becomes weaker until the phase function becomes completely flat
for :math:`g = 0`. For positive values of :math:`g`, the phase function exhibits
an increasingly strong forward scattering peak.


.. figure:: ../../bin/henyey_greenstein.png
   :alt: Henyey-Greenstein phase function

   The Henyey-Greenstein phase function for different values of the asymmetry parameter
   :math:`g`.

Benchmark results
~~~~~~~~~~~~~~~~~

The total absorbed intensity for the different parts of the simulated domain are
given in the table below. The values deviate from the reference benchmark values
given in the exercise. It remains unclear what causes the devaitions but the
stepping scheme remains a possible explanation.

+---------------+----------+-----------------+
| Lower surface |  Medium  |  Upper boundary | 
+===============+==========+=================+
|        0.108  | 0.835    |          0.057  |
+---------------+----------+-----------------+

Part 2
======

For the second part of the exercise, the effect of the parameter :math:`g` on the
radiation field was investigated by repeating the simulation described above for
a range of values of :math:`g` within the interval :math:`[-1, 1]`. The obtained
results are displayed in terms of the fractions of the intensity absorbed in the
different parts of the domain in the figure below.

The intensity leaving the domain through the upper boundary exhibits a peak for
the smallest values of :math:`g`. At the same time, the intensity absorbed by
the medium and the intensity absorbed by the lower surface reach their minima.
As :math:`g` increases, both the intensity absorbed by the medium and that
absorbed by the lower surface increase. However, due to the large optical depth
of the medium, the intensity absorbed by the lower surface remains much smaller
than that absorbed by the medium. At the same time as the intensity absorbed by
the medium and the lower surface decreases, the intensity that leaves the domain
through the upper boundary is reduced.

Also shown in the figure are the simulation results that were obtained with a
Rayleigh-Scattering function. The Rayleigh scattering function is symmetric and
thus has an asymmetry factor of :math:`g = 0`. The results obtained with the
Rayleigh phase function are very close to those obtained with the
Henyey-Greenstein function. This indicates that for :math:`g = 0`, the
Henyey-Greenstein function represents the relevant properties of the Rayleigh
phase function well.


.. figure:: ../../bin/results_4_b.png
   :alt: Simulation results for domain with Henyey-Greenstein function

   Fractions of intensity absorbed in different parts of the simulated domain
   for varying values of the parameter :math:`g` of the Henyey-Greenstein phase
   function. Crosses show the simulation results obtained with the Rayleigh
   scattering function.
