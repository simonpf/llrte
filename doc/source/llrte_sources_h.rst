Radiation sources
=================

The :code:`llrte/sources.h` header provides source classes that can be
used to generate photons in forward or backward Monte Carlo simulations.

Forward sources
---------------

Forward source provide a :code:`sample_photon(generator, atmosphere)` method
that generates a photon and places it on the atmosphere grid.


.. doxygenclass:: llrte::PointSource
   :members:

.. doxygenclass:: llrte::BeamSource
   :members:

Source modifiers
^^^^^^^^^^^^^^^^

In addition to that, a number of source modifiers are provided which manipulate photons after
their creation. These are currently thought for testing purposes only.

.. doxygenclass:: llrte::RandomOffset
   :members:

.. doxygenclass:: llrte::RandomDirection
   :members:

Backward sources
----------------

Backward source provide a :code:`get_intensity(photon)` method
that returns the spectral irradiance for a given direction of an outgoing
photon.

.. doxygenclass:: llrte::PlanarSource
   :members:
