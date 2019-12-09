Exercise 5
----------

For the fifth exercise, the transfer of radiation through a cloud layer was
considered. To run the simulations and generate the figures simply build
the :code:`exercise_5` target using :code:`make`:

.. code:: bash

 $ make exercise_5

Part 1
======

In the first part of the exercise a horizontally-homogeneous cloud layer was
considered. The cloud extends from an altitude of :math:`z_0 = 2` km up to
:math:`z_1 = 7` km. Its optical properties are modeled as

  - Single scattering albedo :math:`a_\text{cloud} = 0.9`
  - Absorption coefficient  :math:`k_\text{cloud} = 1\ \text{km}^{-1}`

For the clear-sky part of the domain the following optical properties are
assumed

  - Single scattering albedo :math:`a_\text{clear} = 1.0`
  - Absorption coefficient  :math:`k_\text{clear} = 2.5\cdot 10^{-3} \text{km}^{-1}`

For the surface an albedo of :math:`A = 0.07` was assumed.

The absorbed intensities in the different parts of the domain are displayed
in the figure below. Since the single-scattering albedo of the clear-sky part
of the domain is 1, no radiation is absorbed outside of the cloud. Within the
cloud, most of the radiation is absorbed at cloud top and decreases towards cloud
base. Since the scene is symmetric along the :math:`x`-axis, no the results
exhibit no horizontal structure except for the MC noise.

.. figure:: ../../bin/results_5_a.png
   :alt: 1D cloud model

   Absorbed intensity for a cloud layer in a horizontally homogeneous
   scene.

The table below displays the fractions of the total intensity that are absorbed
in the different parts of the domain. Due to the low albedo of the surface and
because of the low single-scattering albedo of the cloud, the largest part of
the radiation is absorbed by the surface although it is completely covered by
the cloud. The radiation absorbed by the medium is slightly less than that
absorbed by the surface, while only about :math:`5 \%` of the radiation leave
the domain through the upper boundary.

+---------------+----------+-----------------+
| Lower surface |  Medium  |  Upper boundary | 
+===============+==========+=================+
|        0.516  | 0.439    |          0.046  |
+---------------+----------+-----------------+

Part 2
~~~~~~

For the second part of the exercise, a horizontally heterogeneous cloud
has been simulated. The results are displayed in the figure below. The results
are mostly similar to those obtained in the first part except for the lower flanks
of the cloud, where the absorption is slightly reduced.

.. figure:: ../../bin/results_5_b.png
   :alt: 2D cloud model

   Absorbed intensity for of a horizontally heterogeneous cloud layer. 


The effect of the horizontally heterogeneous cloud layer on the distribution of
the absorbed intensities is displayed in the table below. Since the surface is
only partially covered by cloud, the energy which is absorbed by the surface
increases. At the same time, energy absorbed in the medium decreases compared to
Part 1 since absorption occurs only in the cloud. The radiation that leaves the
domain increases as well, probably because of the strong forward scattering peak
of the Henyey-Greenstein phase function which is used to model the cloud
hydrometeors.

+---------------+----------+-----------------+
| Lower surface |  Medium  |  Upper boundary | 
+===============+==========+=================+
|        0.784  | 0.150    |          0.07   |
+---------------+----------+-----------------+
