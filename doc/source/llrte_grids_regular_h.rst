Grids
=====

An implementation of regular, Cartesian grids is provided in
:code:`llrte/grids/regular.h`. It provides the :code:`RegularGrid` class
which implements the grid and a :code:`GridPosition` class as base class
for objects moving on the grid.

RegularGrid
-----------

A regular grid is described by three arrays, which hold the grid
boundaries in each dimension. It provides the a :code:`GridPosition` struct,
which is a base class for object moving on the grid.

.. doxygenclass:: llrte::RegularGrid
   :members:

GridPosition
------------

The GridPosition struct acts as a base class for object moving on the grid.
It holds their current position and direction as well as the indices of
the grid boundaries that the object will intersect next.

.. note ::
  The indices for the cell boundaries as 1-based. This is used to encode
  leaving particles, which will have index :code:`0` when they leave the
  grid at the lower end.

.. doxygenstruct:: llrte::GridPosition
   :members:
