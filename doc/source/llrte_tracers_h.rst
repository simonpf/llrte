Tracers
=======

Tracer classes provide modular entry points to extract information from Monte
Carlo simulations. Depending on the current event, the Monte Carlo solver calls
a different method of its tracer object. The tracer object is defined by the
NoTrace class which implements all required functions as NOOPs. This is class
is therefore useful as base class for tracers as well as to illustrate the
tracer interface.


NoTrace
-------

As described above, the NoTrace class performs no tracing at all but provides
definitions of all required tracing functions which do nothing.

.. doxygenstruct:: llrte::tracers::NoTrace
   :members:

Histogram
---------

The histogram class traces occurrence frequencies of photons in the different
grids cells of the atmosphere.

.. doxygenclass:: llrte::tracers::Histogram
   :members:

Absorption Tracer
-----------------

The absorption tracer traces absorbed and scattered energy all across the domain
as well as the number of scattering events photons undergo and through which
boundary photons leave the atmosphere.

.. doxygenclass:: llrte::tracers::AbsorptionTracer
   :members:

Photon Tracer
-------------

The photon tracer traces photons that run out of energy or leave the atmosphere
and stores them in a vector.

.. doxygenclass:: llrte::tracers::AbsorptionTracer
   :members:
